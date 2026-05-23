"""
Compare cluster_variants implementations on a real DB.

Algorithms measured
-------------------
greedy_star       Current production algorithm.
                  Sort variants by neighbour count (descending). The highest-degree
                  variant becomes the cluster representative and claims all its
                  neighbours. Claimed variants are skipped in subsequent iterations.
                  BUG: does not check whether a neighbour is already claimed before
                  adding it, so a variant can appear in multiple output clusters.
                  No transitivity: if A overlaps B and B overlaps C but not A, A and C
                  end up in separate clusters.

greedy_star_fixed Same greedy star algorithm with the membership bug fixed.
                  Skips already-claimed neighbours in the inner loop, ensuring each
                  variant appears in exactly one cluster.  No transitivity.

union_find        Union-Find (disjoint-set) on the overlap graph.
                  Every edge (i,j) in the similarity matrix is unioned, giving full
                  transitive closure: A-B and B-C always merge A, B, C into one
                  component even if A and C don't overlap directly.  Representative is
                  the highest-degree member.  Structurally prevents duplicate membership.

union_find_mst    Union-Find + MST-based pruning.
                  Identical Union-Find first pass, then for any component whose posA
                  spread exceeds PRUNE_SPREAD_THRESHOLD, the Minimum Spanning Tree
                  (Kruskal, edge weight = L1 position distance) is computed and its
                  longest edge is cut, splitting the component in two. Recurses until
                  all components are within the spread threshold.

Usage
-----
    python scripts/bench_cluster.py <db_path>
    python scripts/bench_cluster.py <db_path> --skip-existing   # reuse already-run VCFs
    python scripts/bench_cluster.py <db_path> --analyze-only    # no export runs at all

All outputs go to /Users/kselav/Documents/repos/svdb_bench/.
"""

import argparse
import os
import sys
import time
from collections import Counter, defaultdict
from functools import partial
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import svdb.export_module as em
import svdb.database as database
from svdb.export_module import db_header, export
from svdb.ins_similarity import resolve_ins_seq_threshold

BENCH_DIR = "/Users/kselav/Documents/repos/svdb_bench"
PRUNE_SPREAD_THRESHOLD = 2500


# ── cluster implementations ───────────────────────────────────────────────────

def cluster_greedy_star(variant_dictionary, similarity_matrix):
    """Production algorithm — greedy star, highest-degree-first, no transitivity.
    BUG: adds neighbours without checking if already claimed → duplicate membership."""
    cluster_sizes = [[i, len(similarity_matrix[i])] for i in range(len(variant_dictionary))]
    clusters = []
    for i, _ in sorted(cluster_sizes, key=lambda x: x[1], reverse=True):
        if similarity_matrix[i][0] == -1:
            continue
        cluster_dictionary = {}
        for var in similarity_matrix[i]:
            similarity_matrix[var][0] = -1
            cluster_dictionary[var] = variant_dictionary[var]
        clusters.append([variant_dictionary[i], cluster_dictionary])
    return clusters


def cluster_greedy_star_fixed(variant_dictionary, similarity_matrix):
    """Greedy star with the duplicate-membership bug fixed.
    Skips neighbours that are already claimed, so each variant appears in
    exactly one cluster. No transitivity — same clustering semantics as
    greedy_star, just correct."""
    cluster_sizes = [[i, len(similarity_matrix[i])] for i in range(len(variant_dictionary))]
    clusters = []
    for i, _ in sorted(cluster_sizes, key=lambda x: x[1], reverse=True):
        if similarity_matrix[i][0] == -1:
            continue
        cluster_dictionary = {}
        for var in similarity_matrix[i]:
            if similarity_matrix[var][0] != -1:   # skip already-claimed
                similarity_matrix[var][0] = -1
                cluster_dictionary[var] = variant_dictionary[var]
        clusters.append([variant_dictionary[i], cluster_dictionary])
    return clusters


def cluster_union_find(variant_dictionary, similarity_matrix):
    """Union-Find: full transitive closure of the overlap graph."""
    n = len(variant_dictionary)
    parent = list(range(n))
    rank = [0] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px == py:
            return
        if rank[px] < rank[py]:
            px, py = py, px
        parent[py] = px
        if rank[px] == rank[py]:
            rank[px] += 1

    for i in range(n):
        for neighbor in similarity_matrix[i]:
            union(i, int(neighbor))

    components = {}
    for i in range(n):
        root = find(i)
        components.setdefault(root, []).append(i)

    clusters = []
    for members in components.values():
        rep = max(members, key=lambda x: len(similarity_matrix[x]))
        cluster_dict = {m: variant_dictionary[m] for m in members}
        clusters.append([variant_dictionary[rep], cluster_dict])
    return clusters


def _kruskal_mst(members, edges):
    node_idx = {m: k for k, m in enumerate(members)}
    par = list(range(len(members)))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    mst = []
    for w, i, j in sorted(edges):
        pi, pj = find(node_idx[i]), find(node_idx[j])
        if pi != pj:
            par[pi] = pj
            mst.append((w, i, j))
            if len(mst) == len(members) - 1:
                break
    return mst


def _split_component(members, variant_dictionary, similarity_matrix, threshold):
    if len(members) <= 1:
        return [members]
    posA_vals = [variant_dictionary[m]["posA"] for m in members]
    if max(posA_vals) - min(posA_vals) <= threshold:
        return [members]

    member_set = set(members)
    seen = set()
    edges = []
    for i in members:
        for j_raw in similarity_matrix[i]:
            j = int(j_raw)
            if j not in member_set:
                continue
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            seen.add(key)
            vi, vj = variant_dictionary[i], variant_dictionary[j]
            dist = abs(vi["posA"] - vj["posA"]) + abs(vi["posB"] - vj["posB"])
            edges.append((dist, i, j))

    if not edges:
        return [members]
    mst = _kruskal_mst(members, edges)
    if not mst:
        return [members]

    mst_sorted = sorted(mst, key=lambda e: e[0], reverse=True)
    _, cut_i, cut_j = mst_sorted[0]

    adj = defaultdict(set)
    for _, i, j in mst_sorted[1:]:
        adj[i].add(j)
        adj[j].add(i)

    visited = set()
    stack = [cut_i]
    while stack:
        node = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        stack.extend(adj[node] - visited)

    sub_a = list(visited)
    sub_b = [m for m in members if m not in visited]
    return (_split_component(sub_a, variant_dictionary, similarity_matrix, threshold) +
            _split_component(sub_b, variant_dictionary, similarity_matrix, threshold))


def cluster_union_find_mst(variant_dictionary, similarity_matrix,
                           spread_threshold=PRUNE_SPREAD_THRESHOLD):
    """Union-Find + MST pruning: splits components whose posA spread exceeds threshold."""
    n = len(variant_dictionary)
    parent = list(range(n))
    rank = [0] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px == py:
            return
        if rank[px] < rank[py]:
            px, py = py, px
        parent[py] = px
        if rank[px] == rank[py]:
            rank[px] += 1

    for i in range(n):
        for neighbor in similarity_matrix[i]:
            union(i, int(neighbor))

    components = {}
    for i in range(n):
        root = find(i)
        components.setdefault(root, []).append(i)

    final_components = []
    for members in components.values():
        final_components.extend(
            _split_component(members, variant_dictionary, similarity_matrix, spread_threshold)
        )

    clusters = []
    for members in final_components:
        rep = max(members, key=lambda x: len(similarity_matrix[x]))
        cluster_dict = {m: variant_dictionary[m] for m in members}
        clusters.append([variant_dictionary[rep], cluster_dict])
    return clusters


# ── export runner ─────────────────────────────────────────────────────────────

def make_args(db_path, prefix):
    return SimpleNamespace(
        db=db_path,
        prefix=prefix,
        bnd_distance=2500,
        ins_distance=25,
        ins_svlen_ratio=0.90,
        ins_seq_similarity=None,
        data_profile=None,
        no_ins_seq=False,
        overlap=0.8,
        DBSCAN=False,
        epsilon=500,
        min_pts=2,
        memory=False,
        strip_chr=False,
        max_ins_seq_len=500,
        samples="on",
        version="bench",
    )


def run_export(db_path, out_prefix, cluster_fn):
    original = em.cluster_variants
    em.cluster_variants = cluster_fn
    try:
        args = make_args(db_path, out_prefix)
        args.ins_seq_similarity = resolve_ins_seq_threshold(args)
        with database.DB(db_path) as db:
            sample_IDs = db.query_column('SELECT DISTINCT sample FROM SVDB')
        with open(out_prefix + ".vcf", "w") as f:
            f.write(db_header(args) + "\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(
                "\t".join(sample_IDs)))
        t0 = time.perf_counter()
        export(args, sample_IDs)
        return time.perf_counter() - t0
    finally:
        em.cluster_variants = original


# ── analysis helpers ──────────────────────────────────────────────────────────

def count_data_lines(vcf_path):
    with open(vcf_path) as f:
        return sum(1 for ln in f if not ln.startswith("#"))


def count_duplicate_variants(vcf_path):
    seen = Counter()
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            info = line.split("\t")[7]
            for field in info.split(";"):
                if field.startswith("VARIANTS="):
                    for entry in field[9:].split("|"):
                        if entry:
                            seen[entry] += 1
    dups  = sum(1 for c in seen.values() if c > 1)
    extra = sum(c - 1 for c in seen.values() if c > 1)
    return dups, extra


def occ_distribution(vcf_path):
    occ = Counter()
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            info = line.split("\t")[7]
            for field in info.split(";"):
                if field.startswith("OCC="):
                    occ[int(field[4:])] += 1
    return occ


def line_keys(vcf_path):
    keys = set()
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            keys.add((p[0], p[1], p[4]))
    return keys


def stats_by_type(vcf_path):
    lines_by_type = defaultdict(int)
    seen = defaultdict(Counter)
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            info = parts[7]
            svtype = "UNK"
            for field in info.split(";"):
                if field.startswith("SVTYPE="):
                    svtype = field[7:]
                    break
            lines_by_type[svtype] += 1
            for field in info.split(";"):
                if field.startswith("VARIANTS="):
                    for entry in field[9:].split("|"):
                        if entry:
                            seen[svtype][entry] += 1
    result = {}
    for svtype in lines_by_type:
        dups  = sum(1 for c in seen[svtype].values() if c > 1)
        extra = sum(c - 1 for c in seen[svtype].values() if c > 1)
        result[svtype] = (lines_by_type[svtype], dups, extra)
    return result


def investigate_natural_duplicates(db_path):
    """Count (sample, posA, posB) keys that appear in >1 DB row per SV type.

    These are 'natural duplicates' — the same genomic position called twice
    from the same sample with different SVDB idx values.  Any clustering algo
    that produces disjoint components will put them in separate clusters, which
    count_duplicate_variants then flags as duplicate VARIANTS keys.  This function
    measures how many of the flagged duplicates are explained by the input data
    rather than a clustering bug.
    """
    with database.DB(db_path) as db:
        rows = db.query(
            "SELECT var, sample, posA, posB, COUNT(*) as cnt "
            "FROM SVDB "
            "GROUP BY var, sample, posA, posB "
            "HAVING cnt > 1"
        )
    by_type = defaultdict(lambda: [0, 0])  # {var: [n_keys, n_extra]}
    for var, sample, posA, posB, cnt in rows:
        by_type[var][0] += 1
        by_type[var][1] += cnt - 1
    return dict(by_type)


# ── main ─────────────────────────────────────────────────────────────────────

IMPLS = [
    ("greedy_star",
     cluster_greedy_star,
     "greedy star — highest-degree-first, no transitivity [CURRENT DEFAULT, has membership bug]"),
    ("greedy_star_fixed",
     cluster_greedy_star_fixed,
     "greedy star — bug fixed: skips already-claimed neighbours, no transitivity"),
    ("union_find",
     cluster_union_find,
     "Union-Find transitive closure, no size constraint"),
    ("union_find_mst",
     partial(cluster_union_find_mst, spread_threshold=PRUNE_SPREAD_THRESHOLD),
     f"Union-Find + MST pruning (posA spread > {PRUNE_SPREAD_THRESHOLD} bp triggers split)"),
]
IMPL_NAMES = [name for name, _, _ in IMPLS]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("db_path")
    parser.add_argument("--skip-existing", action="store_true",
                        help="skip export for algorithms whose output VCF already exists")
    parser.add_argument("--analyze-only", action="store_true",
                        help="skip all export runs; analyze pre-existing bench VCFs")
    parser.add_argument("--prune-threshold", type=int, default=PRUNE_SPREAD_THRESHOLD)
    args = parser.parse_args()

    impls = list(IMPLS)
    if args.prune_threshold != PRUNE_SPREAD_THRESHOLD:
        impls[-1] = (
            "union_find_mst",
            partial(cluster_union_find_mst, spread_threshold=args.prune_threshold),
            f"Union-Find + MST pruning (posA spread > {args.prune_threshold} bp)",
        )

    print("Algorithms:")
    for name, _, desc in impls:
        print(f"  {name:<22}  {desc}")
    print()

    # ── investigate natural duplicates in the DB ──────────────────────────────
    print("=== DB NATURAL DUPLICATES (same sample+posA+posB, different idx) ===")
    nat_dups = investigate_natural_duplicates(args.db_path)
    if nat_dups:
        for var in sorted(nat_dups):
            n_keys, n_extra = nat_dups[var]
            print(f"  {var:<6}  {n_keys:>6,} duplicate keys  ({n_extra:>6,} extra rows)")
    else:
        print("  none found")
    print()

    timings, vcf_paths = {}, {}

    for name, fn, desc in impls:
        out = os.path.join(BENCH_DIR, f"bench_cluster_{name}")
        vcf_path = out + ".vcf"

        if args.analyze_only or (args.skip_existing and os.path.exists(vcf_path)):
            action = "analyze-only" if args.analyze_only else "reusing existing VCF"
            print(f"=== {name} ({action}) ===")
            vcf_paths[name] = vcf_path
            timings[name] = None
            continue

        print(f"=== {name} ===", flush=True)
        t = run_export(args.db_path, out, fn)
        timings[name] = t
        vcf_paths[name] = vcf_path
        print(f"  elapsed: {t:.1f}s\n", flush=True)

    # ── summary table ─────────────────────────────────────────────────────────
    print("\n=== SUMMARY ===")
    print(f"{'algorithm':<22} {'time':>8}  {'lines':>8}  {'dup_vars':>9}  {'+extra':>7}  OCC (top 5)")
    print("-" * 95)
    for name, _, _ in impls:
        vcf = vcf_paths[name]
        n           = count_data_lines(vcf)
        dups, extra = count_duplicate_variants(vcf)
        t           = timings[name]
        tstr        = f"{t:.1f}s" if t is not None else "  (cached)"
        occ         = occ_distribution(vcf)
        top         = " ".join(f"{k}×{v}" for k, v in sorted(occ.items(), reverse=True)[:5])
        print(f"{name:<22} {tstr:>10}  {n:>8,}  {dups:>9,}  {extra:>7,}  {top}")

    # ── per-type breakdown ────────────────────────────────────────────────────
    all_types = sorted(set(
        t for name, _, _ in impls for t in stats_by_type(vcf_paths[name])
    ))
    print(f"\n=== PER-TYPE BREAKDOWN ===")
    print(f"{'type':<8} {'algorithm':<22} {'lines':>8}  {'dup_vars':>9}  {'+extra':>7}")
    print("-" * 60)
    for svtype in all_types:
        for name, _, _ in impls:
            s = stats_by_type(vcf_paths[name])
            lines, dups, extra = s.get(svtype, (0, 0, 0))
            print(f"{svtype:<8} {name:<22} {lines:>8,}  {dups:>9,}  {extra:>7,}")
        print()

    # ── diffs vs greedy_star ──────────────────────────────────────────────────
    sq_keys = line_keys(vcf_paths["greedy_star"])
    for name, _, _ in impls:
        if name == "greedy_star":
            continue
        other_keys = line_keys(vcf_paths[name])
        only_sq    = sq_keys - other_keys
        only_other = other_keys - sq_keys
        print(f"\n=== DIFF {name} vs greedy_star ===")
        print(f"  lines only in greedy_star  : {len(only_sq):,}")
        print(f"  lines only in {name:<18}: {len(only_other):,}")
        if only_sq:
            print(f"  sample from greedy_star only (CHROM, POS, ALT):")
            for k in sorted(only_sq)[:3]:
                print(f"    {k}")
        if only_other:
            print(f"  sample from {name} only (CHROM, POS, ALT):")
            for k in sorted(only_other)[:3]:
                print(f"    {k}")


if __name__ == "__main__":
    main()
