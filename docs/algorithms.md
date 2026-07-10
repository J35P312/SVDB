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
   `--ins_seq_similarity` for insertions with known sequence. Sequences over
   10% ambiguous 'N' bases (`ins_similarity.has_excess_n`) are treated as
   absent — N content makes Levenshtein similarity meaningless, so such
   insertions fall back to position+SVLEN matching instead.

The similarity graph (which variants can be merged with which) is then
collapsed by one of two algorithms:

### Greedy star (`--cluster_method star`, default)

```text
Sort variants by degree (overlapping neighbours) descending.
For each variant, highest-degree first:
    if already marked as non-representative: skip
    otherwise: become representative
        for each neighbour:
            mark neighbour as non-representative
            add neighbour to this cluster
```

Marking a variant as non-representative only prevents it from *starting* a
new cluster — it does not exclude it from being collected by other active
representatives whose neighbour list includes it.  A variant that overlaps
two separate representatives therefore contributes to both clusters.

**Properties**: fast; no transitivity (A-B and B-C do not force A and C into
the same cluster). A variant in an ambiguous region between two groups
appears in each — OCC and FRQ reflect every group it genuinely overlaps
rather than being forced into one by an arbitrary assignment.
On a 100-sample SV dataset ~7 % of variants appear in more than one cluster;
the maximum observed is 8 clusters (empirically bounded by local density).

**Use when**: you want OCC to accurately represent all overlapping groups —
typical for population frequency databases.  Use `union_find` if you need
hard exclusive membership (e.g. downstream tools that assume disjoint clusters).

### Union-Find (`--cluster_method union_find`)

Builds the full transitive closure of the overlap graph using a path-compressed
Union-Find (disjoint-set) structure:

```text
For every overlapping pair (i, j):
    union(i, j)
Collect connected components; representative = highest-degree member.
```

**Properties**: A-B and B-C always merge {A, B, C} even if A and C do not
overlap directly. Each variant belongs to exactly one cluster.
Produces fewer, larger clusters with higher OCC counts.

**Use when**: you need hard exclusive membership — e.g. downstream tools that
assume disjoint clusters, or when you prefer larger merged groups over
accurate per-group frequencies.

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
