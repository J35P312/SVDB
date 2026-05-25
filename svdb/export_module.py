import logging
import multiprocessing
import os
import sys
from collections import Counter

import numpy as np

from . import dbscan, database, ins_similarity
from .ins_similarity import decompress_ins_seq
from .vcf_utils import normalize_chrom

logger = logging.getLogger(__name__)

try:
    from scipy.spatial import cKDTree as _cKDTree  # type: ignore[import-untyped]
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False

# Switch to cKDTree above this cluster size; numpy SIMD wins below it.
_KDTREE_MIN_SIZE = 200


def make_representing_variant(variant_type, chrA, chrB, posA, ci_A_start, ci_A_end, posB, ci_B_start, ci_B_end, ins_seq=None, ins_len=None):
    """Build the representing-variant dict used as cluster[0] in vcf_line."""
    return {
        "type": variant_type,
        "chrA": chrA,
        "chrB": chrB,
        "posA": posA,
        "ci_A_start": ci_A_start,
        "ci_A_end": ci_A_end,
        "posB": posB,
        "ci_B_start": ci_B_start,
        "ci_B_end": ci_B_end,
        "ins_seq": ins_seq,
        "ins_len": ins_len,
    }


def build_genotype_columns(sample_IDs, hit_sample_ids):
    """Return per-sample GT strings in sample_IDs order.

    Samples present in hit_sample_ids get './1'; all others get '0/0'.
    """
    zygosity = {sample: "0/0" for sample in sample_IDs}
    for sid in hit_sample_ids:
        zygosity[sid] = "./1"
    return [zygosity[sample] for sample in sample_IDs]


def fetch_all_variants(variant, chrA, chrB, db, max_ins_seq_len=None):
    """Load all variant data for one (var, chrA, chrB) group in a single query.

    Returns (variant_dict, pos_coords, pos_indices):
      variant_dict   – {idx: {posA, posB, sample_id, ins_seq, ins_len}}
      pos_coords     – np.array shape (n, 2) columns [posA, posB]
      pos_indices    – np.array shape (n,) of idx values, parallel to pos_coords
    Returns ({}, empty, empty) when the group has no rows.

    max_ins_seq_len: if set, sequences longer than this are treated as None
    during comparison (position+SVLEN only); the full sequence is still written
    to the ALT field of the output VCF.
    """
    has_ins = db.has_ins_table()
    if has_ins:
        hits = db.query(
            'SELECT s.posA, s.posB, s.sample, s.idx, i.ins_seq, i.ins_len '
            'FROM SVDB s LEFT JOIN INS i ON s.idx = i.idx '
            f'WHERE s.var == \'{variant}\' AND s.chrA == \'{chrA}\' AND s.chrB == \'{chrB}\''
        )
    else:
        hits = db.query(
            f'SELECT posA, posB, sample, idx FROM SVDB '
            f'WHERE var == \'{variant}\' AND chrA == \'{chrA}\' AND chrB == \'{chrB}\''
        )
    if not hits:
        return {}, np.empty((0, 2), dtype=int), np.array([], dtype=int)

    variant_dict = {}
    x, y, indices = [], [], []
    for hit in hits:
        idx = int(hit[3])
        seq = decompress_ins_seq(hit[4]) if has_ins else None
        if seq and max_ins_seq_len is not None and len(seq) > max_ins_seq_len:
            seq = None
        variant_dict[idx] = {
            "posA": int(hit[0]),
            "posB": int(hit[1]),
            "sample_id": hit[2],
            "ins_seq": seq,
            "ins_len": hit[5] if has_ins else None,
        }
        x.append(int(hit[0]))
        y.append(int(hit[1]))
        indices.append(idx)
    return variant_dict, np.column_stack((x, y)), np.array(indices)


def db_header(args):
    no_samples = getattr(args, "samples", "on") == "off"
    headerString = '##fileformat=VCFv4.1\n'
    headerString += '##source=SVDB\n'
    headerString += '##ALT=<ID=DEL,Description="Deletion">\n'
    headerString += '##ALT=<ID=DUP,Description="Duplication">\n'
    headerString += '##ALT=<ID=INV,Description="Inversion">\n'
    headerString += '##ALT=<ID=INS,Description="Insertion">\n'
    headerString += '##ALT=<ID=BND,Description="Break end">\n'
    headerString += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
    headerString += '##INFO=<ID=END,Number=1,Type=String,Description="End of an intra-chromosomal variant">\n'
    headerString += '##INFO=<ID=OCC,Number=1,Type=Integer,Description="The number of occurences of the event in the database">\n'
    headerString += '##INFO=<ID=NSAMPLES,Number=1,Type=Integer,Description="the number of samples within the database">\n'
    headerString += '##INFO=<ID=VARIANTS,Number=1,Type=Integer,Description="a| separated list of the positions of the clustered variants">\n'
    headerString += '##INFO=<ID=FRQ,Number=1,Type=Float,Description="the frequency of the variant">\n'
    headerString += '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
    headerString += '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n'
    headerString += '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">\n'
    if not no_samples:
        headerString += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    headerString += '##SVDB_version={} cmd=\"{}\"'.format(args.version, " ".join(sys.argv))
    return headerString


def vcf_line(cluster, id_tag, sample_IDs, strip_chr=False, no_samples=False):
    _chrom = normalize_chrom if strip_chr else (lambda x: x)
    info_field = "SVTYPE={};".format(cluster[0]["type"])
    vcf_line = []
    vcf_line.append(_chrom(cluster[0]["chrA"]))
    vcf_line.append(str(cluster[0]["posA"]))
    vcf_line.append(id_tag)
    vcf_line.append("N")
    is_ins = cluster[0]["type"] == "INS" and cluster[0]["chrA"] == cluster[0]["chrB"]
    if cluster[0]["chrA"] == cluster[0]["chrB"] and cluster[0]["type"] != "BND":
        ins_seq = cluster[0].get("ins_seq")
        ins_len = cluster[0].get("ins_len")
        if is_ins and ins_seq:
            vcf_line.append("N" + ins_seq)
            info_field += "SVLEN={};".format(len(ins_seq))
        elif is_ins and ins_len:
            vcf_line.append("<INS>")
            info_field += "END={};SVLEN={};".format(cluster[0]["posB"], ins_len)
        else:
            vcf_line.append("<" + cluster[0]["type"] + ">")
            info_field += "END={};SVLEN={};".format(cluster[0]["posB"], abs(cluster[0]["posA"] - cluster[0]["posB"]))
    else:
        vcf_line.append("N[{}:{}[".format(_chrom(cluster[0]["chrB"]), cluster[0]["posB"]))

    sample_set = set([])
    CIPOS = []
    CIEND = []
    for variant in cluster[1]:
        CIPOS.append(cluster[1][variant]["posA"])
        CIEND.append(cluster[1][variant]["posB"])
        sample_set.add(cluster[1][variant]["sample_id"])

    CIPOS_start = -abs(cluster[0]["posA"] - min(CIPOS))
    CIPOS_end = abs(cluster[0]["posA"] - max(CIPOS))

    CIEND_start = -abs(cluster[0]["posB"] - min(CIEND))
    CIEND_end = abs(cluster[0]["posB"] - max(CIEND))

    info_field += "NSAMPLES={};OCC={};FRQ={};CIPOS={},{};CIEND={},{};".format(len(sample_IDs), len(
        sample_set), round(len(sample_set) / float(len(sample_IDs)), 4), CIPOS_start, CIPOS_end, CIEND_start, CIEND_end)
    variant_field = "VARIANTS="
    for variant in cluster[1]:
        if no_samples:
            variant_field += "|{}:{}".format(cluster[1][variant]["posA"], cluster[1][variant]["posB"])
        else:
            variant_field += "|{}:{}:{}".format(cluster[1][variant]["sample_id"], cluster[1][variant]["posA"], cluster[1][variant]["posB"])
    info_field += variant_field
    vcf_line.append(".")
    vcf_line.append("PASS")
    vcf_line.append(info_field)
    if not no_samples:
        hit_sample_ids = [cluster[1][v]["sample_id"] for v in cluster[1]]
        format_cols = build_genotype_columns(sample_IDs, hit_sample_ids)
        vcf_line.append("GT")
        vcf_line.append("\t".join(format_cols))
    return "\t".join(vcf_line)


def expand_chain(chain, coordinates, chrA, chrB, distance, overlap,
                 ins_svlen_ratio=None, ins_seq_threshold=None, no_ins_seq=False):
    """Build pairwise match lists for every variant in the cluster.

    `chain` is a dict {0..n-1: variant_dict} and `coordinates` is an (n,3)
    int array with columns [chain_idx, posA, posB].  Because both are built
    with the same enumerate(), coordinates[k,0]==k for all k, so row indices
    and chain keys are interchangeable — this is exploited throughout.

    Overlap gates are evaluated in cheapest-first order:
      1. Spatial lookup    — numpy mask or cKDTree (L∞), already vectorised
      2. Overlap ratio     — vectorised numpy (DEL/DUP/INV only)
      3. SVLEN ratio       — vectorised numpy (INS only)
      4. Sequence gate     — per-pair rapidfuzz (INS only, irreducibly serial)
    """
    is_ins = chrA == chrB and overlap == -1
    n = len(chain)
    chain_data = {}

    use_tree = _SCIPY_AVAILABLE and n >= _KDTREE_MIN_SIZE
    tree = _cKDTree(coordinates[:, 1:3].astype(float)) if use_tree else None

    # Precompute per-variant arrays so the inner loop avoids repeated dict access.
    # Safe because chain keys are exactly 0..n-1.
    if is_ins and ins_svlen_ratio is not None:
        all_lens = np.array([chain[k].get("ins_len") or 0 for k in range(n)],
                            dtype=np.float64)

    need_seq = is_ins and not no_ins_seq and ins_seq_threshold is not None
    if need_seq:
        all_seqs = [chain[k].get("ins_seq") or "" for k in range(n)]

    for i in range(n):
        variant = chain[i]
        va_posA = variant["posA"]
        va_posB = variant["posB"]

        # ── 1. spatial lookup ────────────────────────────────────────────────
        # Returns row indices into `coordinates`, which equal chain keys.
        # L∞ norm: max(|dposA|, |dposB|) ≤ distance — identical to both
        # precise_overlap (INS/BND) and the position guard in isSameVariation.
        if use_tree:
            cand_idxs = np.array(
                tree.query_ball_point([va_posA, va_posB], distance, p=np.inf),
                dtype=np.int64)
        else:
            cand_idxs = np.where(
                (np.abs(coordinates[:, 1] - va_posA) <= distance) &
                (np.abs(coordinates[:, 2] - va_posB) <= distance)
            )[0].astype(np.int64)

        if len(cand_idxs) == 0:
            chain_data[i] = np.array([], dtype=np.int64)
            continue

        # ── 2. overlap ratio (DEL / DUP / INV) ──────────────────────────────
        # INS and BND: position check already equals the full gate → all pass.
        if chrA != chrB or is_ins:
            match_mask = np.ones(len(cand_idxs), dtype=bool)
        else:
            cand_posA = coordinates[cand_idxs, 1]
            cand_posB = coordinates[cand_idxs, 2]
            region_start  = np.minimum(cand_posA, va_posA)
            overlap_start = np.maximum(cand_posA, va_posA)
            region_end    = np.maximum(cand_posB, va_posB)
            overlap_end   = np.minimum(cand_posB, va_posB)
            denom = (region_end - region_start + 1).astype(np.float64)
            ratio_vals = (overlap_end - overlap_start + 1) / denom
            match_mask = ratio_vals >= overlap

        # ── 3. SVLEN ratio (INS, vectorised) ────────────────────────────────
        if is_ins and ins_svlen_ratio is not None:
            va_len = variant.get("ins_len") or 0
            if va_len > 0:
                cand_lens = all_lens[cand_idxs]
                valid = cand_lens > 0
                svlen_ratio = np.where(
                    valid,
                    np.minimum(cand_lens, va_len) / np.maximum(cand_lens, va_len),
                    1.0,
                )
                match_mask &= svlen_ratio >= ins_svlen_ratio

        # ── 4. sequence gate (INS, per-pair) ────────────────────────────────
        if need_seq:
            va_seq = variant.get("ins_seq") or ""
            for k in np.where(match_mask)[0]:
                if not ins_similarity.sequence_gate(
                        va_seq, all_seqs[int(cand_idxs[k])], ins_seq_threshold):
                    match_mask[k] = False

        chain_data[i] = cand_idxs[match_mask]
    return chain_data


def cluster_variants(variant_dictionary, similarity_matrix):
    """Greedy star clustering (default).

    Highest-degree variant becomes representative and claims its neighbours.
    A variant that overlaps multiple representatives contributes to each of
    their clusters — intentional, so OCC/FRQ reflect all groups it belongs to.
    No transitivity: A-B and B-C do not force A and C into the same cluster.
    Use --cluster_method union_find for exclusive (hard) cluster membership.
    """
    cluster_sizes = [[i, len(similarity_matrix[i])] for i in range(len(variant_dictionary))]

    clusters = []
    for i, _ in sorted(cluster_sizes, key=lambda x: (x[1]), reverse=True):
        if similarity_matrix[i][0] == -1:
            continue

        cluster_dictionary = {}
        for var in similarity_matrix[i]:
            similarity_matrix[var][0] = -1
            cluster_dictionary[var] = variant_dictionary[var]
        variant = variant_dictionary[i]

        clusters.append([variant, cluster_dictionary])
    return clusters


def cluster_variants_union_find(variant_dictionary, similarity_matrix):
    """Union-Find clustering.

    Full transitive closure of the overlap graph: A-B and B-C always merge
    A, B, C into one component even if A and C don't overlap directly.
    Structurally prevents duplicate membership.  Representative is the
    highest-degree member of each component.
    """
    n = len(variant_dictionary)
    parent = list(range(n))
    rank = [0] * n

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
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

    components: dict[int, list[int]] = {}
    for i in range(n):
        components.setdefault(find(i), []).append(i)

    clusters = []
    for members in components.values():
        rep = max(members, key=lambda x: len(similarity_matrix[x]))
        cluster_dict = {m: variant_dictionary[m] for m in members}
        clusters.append([variant_dictionary[rep], cluster_dict])
    return clusters


def _pick_ins_seq(variant_dict):
    """Return the most common non-null ins_seq across the cluster, or None."""
    seqs = [v.get("ins_seq") for v in variant_dict.values() if v.get("ins_seq")]
    if not seqs:
        return None
    return Counter(seqs).most_common(1)[0][0]


def _pick_ins_len(variant_dict):
    """Return the most common non-null ins_len across the cluster, or None."""
    lens = [v.get("ins_len") for v in variant_dict.values() if v.get("ins_len")]
    if not lens:
        return None
    return Counter(lens).most_common(1)[0][0]


def _expand_and_cluster_worker(task: tuple) -> list:
    """Worker: run expand_chain + cluster for one spatial cluster group.

    All arguments are plain scalars, dicts, or ndarrays so the tuple
    pickles cleanly for multiprocessing.Pool dispatch.
    """
    (sub_dict, sub_coords, variant, chrA, chrB,
     ins_distance, bnd_distance, overlap,
     ins_svlen_ratio, ins_seq_threshold, no_ins_seq, cluster_method) = task

    if "INS" in variant:
        similarity_matrix = expand_chain(
            sub_dict, sub_coords, chrA, chrB, ins_distance, -1,
            ins_svlen_ratio=ins_svlen_ratio,
            ins_seq_threshold=ins_seq_threshold,
            no_ins_seq=no_ins_seq,
        )
    else:
        similarity_matrix = expand_chain(
            sub_dict, sub_coords, chrA, chrB, bnd_distance, overlap)

    cluster_fn = (cluster_variants_union_find
                  if cluster_method == "union_find"
                  else cluster_variants)
    clusters = cluster_fn(sub_dict, similarity_matrix)

    result = []
    for rep, cdict in clusters:
        rep["type"] = variant
        rep["chrA"] = chrA
        rep["chrB"] = chrB
        rep["ins_seq"] = _pick_ins_seq(cdict)
        rep["ins_len"] = _pick_ins_len(cdict)
        result.append((rep, cdict))
    return result


def _resolve_workers(requested: int) -> int:
    """0 = auto (all logical CPUs); ≥1 = explicit count."""
    if requested <= 0:
        return max(1, os.cpu_count() or 1)
    return requested


def _run_tasks(tasks: list, workers: int) -> list:
    """Dispatch tasks to a process pool; fall back to serial on any failure."""
    if not tasks:
        return []
    w = min(workers, len(tasks))
    if w <= 1:
        return [_expand_and_cluster_worker(t) for t in tasks]
    logger.debug("clustering %d group(s) with %d worker(s)", len(tasks), w)
    chunksize = max(1, len(tasks) // (w * 4))
    try:
        ctx = multiprocessing.get_context("spawn")
        with ctx.Pool(processes=w) as pool:
            return pool.map(_expand_and_cluster_worker, tasks, chunksize=chunksize)
    except Exception as exc:
        logger.warning("parallel clustering failed (%s); using serial fallback", exc)
        return [_expand_and_cluster_worker(t) for t in tasks]


def svdb_cluster_main(chrA, chrB, variant, sample_IDs, args, db, i, f):
    workers = _resolve_workers(getattr(args, "workers", 1))
    max_ins_seq_len = getattr(args, "max_ins_seq_len", None)
    all_data, pos_coords, pos_indices = fetch_all_variants(variant, chrA, chrB, db, max_ins_seq_len)
    if not all_data:
        return i

    if args.coarse:
        labels = dbscan.main(pos_coords, args.epsilon, args.min_pts)
    elif "INS" in variant:
        labels = dbscan.main(pos_coords, args.ins_distance, 2)
    else:
        labels = dbscan.main(pos_coords, args.bnd_distance, 2)

    strip_chr = getattr(args, "strip_chr", False)
    no_samples = getattr(args, "samples", "on") == "off"
    unique_labels = set(labels)

    # singletons: written immediately, no clustering needed
    singleton_mask = labels == -1
    for pos, idx in zip(pos_coords[singleton_mask], pos_indices[singleton_mask]):
        v = all_data[idx]
        representing_var = make_representing_variant(
            variant, chrA, chrB, pos[0], pos[0], pos[0], pos[1], pos[1], pos[1],
            v.get("ins_seq"), v.get("ins_len"))
        f.write(vcf_line([representing_var, {0: v}], f"cluster_{i}", sample_IDs, strip_chr, no_samples) + "\n")
        i += 1

    # non-singleton clusters
    tasks = []
    for unique_label in unique_labels:
        if unique_label == -1:
            continue
        mask = labels == unique_label
        cluster_pos = pos_coords[mask]
        cluster_idxs = pos_indices[mask]
        sub_dict = {j: all_data[idx] for j, idx in enumerate(cluster_idxs)}
        sub_coords = np.array([[j, all_data[idx]["posA"], all_data[idx]["posB"]]
                                for j, idx in enumerate(cluster_idxs)])

        if args.coarse:
            avg_point = np.array([np.mean(cluster_pos[:, 0]), np.mean(cluster_pos[:, 1])])
            ins_seq = _pick_ins_seq(sub_dict)
            ins_len = _pick_ins_len(sub_dict)
            representing_var = make_representing_variant(
                variant, chrA, chrB,
                int(avg_point[0]), np.amin(cluster_pos[:, 0]), np.amax(cluster_pos[:, 0]),
                int(avg_point[1]), np.amin(cluster_pos[:, 1]), np.amax(cluster_pos[:, 1]),
                ins_seq, ins_len)
            f.write(vcf_line([representing_var, sub_dict], f"cluster_{i}", sample_IDs, strip_chr, no_samples) + "\n")
            i += 1
        else:
            tasks.append((
                sub_dict, sub_coords, variant, chrA, chrB,
                args.ins_distance, args.bnd_distance, args.overlap,
                args.ins_svlen_ratio,
                args.ins_seq_similarity,
                args.no_ins_seq,
                args.cluster_method,
            ))

    for group_clusters in _run_tasks(tasks, workers):
        for rep, cdict in group_clusters:
            f.write(vcf_line([rep, cdict], f"cluster_{i}", sample_IDs, strip_chr, no_samples) + "\n")
            i += 1

    return i


def export(args, sample_IDs):
    ins_similarity.apply_ins_profile(args)
    with database.DB(args.db, memory=args.memory) as db:
        groups = db.query('SELECT DISTINCT var, chrA, chrB FROM SVDB')

        if any("INS" in var for var, _, _ in groups) and not db.has_ins_table():
            logger.warning(
                "database does not contain insertion sequence/length data — "
                "exporting insertions without sequence in ALT column. "
                "To enable full insertion export, run: svdb --build --upgrade --files <original_vcfs>"
            )
        elif (any("INS" in var for var, _, _ in groups)
              and db.has_ins_table()
              and getattr(args, "max_ins_seq_len", None) is None
              and not getattr(args, "no_ins_seq", False)):
            logger.info(
                "No --max_ins_seq_len set; sequence similarity will run on all insertion lengths. "
                "For large databases, --max_ins_seq_len 500 or --max_ins_seq_len 1000 "
                "is recommended to significantly reduce export time."
            )

        i = 0
        with open(args.prefix + ".vcf", 'a') as f:
            for variant, chrA, chrB in groups:
                i = svdb_cluster_main(chrA, chrB, variant, sample_IDs, args, db, i, f)


def main(args):
    sample_IDs = []
    if not args.prefix:
        args.prefix = args.db.replace(".db", "")

    with database.DB(args.db) as db:
        sample_IDs = db.query_column('SELECT DISTINCT sample FROM SVDB')

    no_samples = getattr(args, "samples", "on") == "off"
    with open(args.prefix + ".vcf", 'w') as f:
        f.write(db_header(args) + "\n")
        if no_samples:
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        else:
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(sample_IDs)))
    export(args, sample_IDs)
