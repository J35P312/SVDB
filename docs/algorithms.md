# SVDB export clustering algorithms

SVDB export uses a two-pass pipeline to cluster database variants into representative VCF records.

## Overview

```text
All variants (chrA, chrB, type)
    └─ Pass 1: DBSCAN spatial grouping
           groups nearby variants into candidate clusters
           ↓
       singletons (noise) → written immediately
       ↓
    └─ Pass 2: second-pass refinement   [skipped with --coarse]
           overlap ratio + SVLEN ratio + sequence similarity gates
           → similarity graph
           → star or union_find clustering
           → one VCF line per cluster
```

---

## Pass 1: DBSCAN spatial grouping

**What it does**: groups variants by genomic proximity using a sliding-window
density approach inspired by DBSCAN (Density-Based Spatial Clustering of
Applications with Noise; Ester et al., KDD 1996).

**Parameters**:

| Flag | Role | Default |
| ---- | ---- | ------- |
| `--epsilon` | spatial radius (bp): variants within this distance seed a cluster | 500 |
| `--min_pts` | minimum consecutive variants within epsilon to form a cluster | 2 |

For insertions the radius is `--ins_distance` (default 25, or 50 for
`--data_profile cohort`/`position_only`). For other SV types it is
`--bnd_distance`.

**Variants labeled as noise** (isolated, no neighbors within epsilon) become
singletons and are written directly to the output without entering Pass 2.

---

## Pass 2: second-pass refinement (default)

Within each DBSCAN spatial group, pairwise overlap is evaluated through a
cascade of gates — cheapest first to minimise work:

1. **Spatial gate** (L∞ distance, vectorised): enforces `--ins_distance` /
   `--bnd_distance` / `--overlap` as a hard window.
2. **SVLEN ratio gate** (vectorised): requires `min(svlen)/max(svlen) ≥ --ins_svlen_ratio`
   for insertions with known length.
3. **Sequence similarity gate** (Levenshtein, per-pair): requires similarity ≥
   `--ins_seq_similarity` for insertions with known sequence.

The similarity graph (which variants can be merged with which) is then
collapsed by one of two algorithms:

### Greedy star (`--cluster_method star`, default)

```text
Sort variants by degree (overlapping neighbours) descending.
While unclaimed variants remain:
    pick the highest-degree unclaimed variant as representative
    claim all its unclaimed neighbours
```

**Properties**: fast; no transitivity. If A overlaps B and B overlaps C but A
does not overlap C, star keeps them in separate clusters. Produces more,
smaller clusters.

**Use when**: you want conservative merging — every cluster member directly
overlaps the representative.

### Union-Find (`--cluster_method union_find`)

Builds the full transitive closure of the overlap graph using a path-compressed
Union-Find (disjoint-set) structure:

```text
For every overlapping pair (i, j):
    union(i, j)
Collect connected components; representative = highest-degree member.
```

**Properties**: A-B and B-C always merge {A, B, C} even if A and C do not
overlap directly. Prevents duplicate cluster membership by construction.
Produces fewer, larger clusters with higher OCC counts.

**Use when**: maximal merging across a cohort; population-level databases where
transitive relationships are expected.

---

## `--coarse` mode

Skips Pass 2 entirely. Uses DBSCAN spatial groups directly:
representative = centroid of the group; ALT = most-common insertion sequence
across members. Faster but less precise; controlled by `--epsilon` and
`--min_pts` only.

`--DBSCAN` is a deprecated alias for `--coarse` (accepted, emits a warning).

---

## Parallelism

Pass 2 is parallelised: each DBSCAN group is dispatched as an independent task
to a `multiprocessing.Pool` (spawn context, safe with numpy on all platforms).
Workers receive plain dicts and arrays — no file handles or DB connections.
Results are written in deterministic order.

`--workers 0` (default): use all logical CPUs.
`--workers 1`: serial; useful for profiling or shared-system courtesy.

See [architecture.md](architecture.md) for implementation-level detail.
