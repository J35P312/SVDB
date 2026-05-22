import logging
import sys
from collections import Counter

import numpy as np

from . import dbscan, database, ins_similarity, overlap_module
from .vcf_utils import normalize_chrom

logger = logging.getLogger(__name__)


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


def fetch_all_variants(variant, chrA, chrB, db):
    """Load all variant data for one (var, chrA, chrB) group in a single query.

    Returns (variant_dict, pos_coords, pos_indices):
      variant_dict   – {idx: {posA, posB, sample_id, ins_seq, ins_len}}
      pos_coords     – np.array shape (n, 2) columns [posA, posB]
      pos_indices    – np.array shape (n,) of idx values, parallel to pos_coords
    Returns ({}, empty, empty) when the group has no rows.
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
        variant_dict[idx] = {
            "posA": int(hit[0]),
            "posB": int(hit[1]),
            "sample_id": hit[2],
            "ins_seq": hit[4] if has_ins else None,
            "ins_len": hit[5] if has_ins else None,
        }
        x.append(int(hit[0]))
        y.append(int(hit[1]))
        indices.append(idx)
    return variant_dict, np.column_stack((x, y)), np.array(indices)


def db_header(args):
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
    headerString += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    headerString += '##SVDB_version={} cmd=\"{}\"'.format(args.version, " ".join(sys.argv))
    return headerString


def vcf_line(cluster, id_tag, sample_IDs, strip_chr=False):
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
        variant_field += "|{}:{}:{}".format(cluster[1][variant]["sample_id"], cluster[1][variant]["posA"], cluster[1][variant]["posB"])
    info_field += variant_field
    vcf_line.append(".")
    vcf_line.append("PASS")
    vcf_line.append(info_field)
    hit_sample_ids = [cluster[1][v]["sample_id"] for v in cluster[1]]
    format_cols = build_genotype_columns(sample_IDs, hit_sample_ids)
    vcf_line.append("GT")
    vcf_line.append("\t".join(format_cols))
    return "\t".join(vcf_line)


def expand_chain(chain, coordinates, chrA, chrB, distance, overlap,
                 ins_svlen_ratio=None, ins_seq_threshold=None, no_ins_seq=False):
    is_ins = chrA == chrB and overlap == -1
    chain_data = {}
    for i, idx in enumerate(chain):
        chain_data[i] = []
        variant = chain[idx]

        rows = coordinates[(distance >= abs(coordinates[:, 1] - variant["posA"]))
                           & (distance >= abs(coordinates[:, 2] - variant["posB"]))]

        candidates = rows[:, 0]
        for candidate in candidates:
            var = chain[candidate]
            if chrA != chrB:
                match = True
            elif is_ins:
                _, match = overlap_module.precise_overlap(
                    variant["posA"], variant["posB"], var["posA"], var["posB"], distance)
            else:
                _, match = overlap_module.isSameVariation(
                    variant["posA"], variant["posB"], var["posA"], var["posB"], overlap, distance)

            if match and is_ins and ins_svlen_ratio is not None:
                len_a = variant.get("ins_len")
                len_b = var.get("ins_len")
                if len_a is not None and len_b is not None:
                    if not overlap_module.insertion_svlen_match(len_a, len_b, ins_svlen_ratio):
                        match = False

            if match and is_ins and not no_ins_seq and ins_seq_threshold is not None:
                seq_a = variant.get("ins_seq") or ""
                seq_b = var.get("ins_seq") or ""
                if not ins_similarity.sequence_gate(seq_a, seq_b, ins_seq_threshold):
                    match = False

            if match:
                chain_data[i].append(candidate)

        chain_data[i] = np.array(chain_data[i])
    return chain_data


def cluster_variants(variant_dictionary, similarity_matrix):
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


def overlap_cluster(variant_dictionary, coordinates, variant, chrA, chrB, sample_IDs, args, f, i):
    if "INS" in variant:
        similarity_matrix = expand_chain(
            variant_dictionary, coordinates, chrA, chrB, args.ins_distance, -1,
            ins_svlen_ratio=getattr(args, "ins_svlen_ratio", None),
            ins_seq_threshold=getattr(args, "ins_seq_similarity", None),
            no_ins_seq=getattr(args, "no_ins_seq", False))
    else:
        similarity_matrix = expand_chain(
           variant_dictionary, coordinates, chrA, chrB, args.bnd_distance, args.overlap)

    strip_chr = getattr(args, "strip_chr", False)
    clusters = cluster_variants(variant_dictionary, similarity_matrix)
    for clustered_variants in clusters:
        clustered_variants[0]["type"] = variant
        clustered_variants[0]["chrA"] = chrA
        clustered_variants[0]["chrB"] = chrB
        clustered_variants[0]["ins_seq"] = _pick_ins_seq(clustered_variants[1])
        clustered_variants[0]["ins_len"] = _pick_ins_len(clustered_variants[1])
        f.write(vcf_line(clustered_variants, f"cluster_{i}", sample_IDs, strip_chr) + "\n")
        i += 1
    return i


def svdb_cluster_main(chrA, chrB, variant, sample_IDs, args, db, i, f):
    all_data, pos_coords, pos_indices = fetch_all_variants(variant, chrA, chrB, db)
    if not all_data:
        return i

    if args.DBSCAN:
        labels = dbscan.main(pos_coords, args.epsilon, args.min_pts)
    elif "INS" in variant:
        labels = dbscan.main(pos_coords, args.ins_distance, 2)
    else:
        labels = dbscan.main(pos_coords, args.bnd_distance, 2)

    strip_chr = getattr(args, "strip_chr", False)
    unique_labels = set(labels)

    # singletons
    singleton_mask = labels == -1
    for pos, idx in zip(pos_coords[singleton_mask], pos_indices[singleton_mask]):
        v = all_data[idx]
        variant_dictionary = {0: v}
        representing_var = make_representing_variant(
            variant, chrA, chrB, pos[0], pos[0], pos[0], pos[1], pos[1], pos[1],
            v.get("ins_seq"), v.get("ins_len"))
        cluster = [representing_var, variant_dictionary]
        f.write(vcf_line(cluster, f"cluster_{i}", sample_IDs, strip_chr) + "\n")
        i += 1

    # clusters
    for unique_label in unique_labels:
        if unique_label == -1:
            continue
        mask = labels == unique_label
        cluster_pos = pos_coords[mask]
        cluster_idxs = pos_indices[mask]

        sub_dict = {j: all_data[idx] for j, idx in enumerate(cluster_idxs)}
        sub_coords = np.array([[j, all_data[idx]["posA"], all_data[idx]["posB"]]
                                for j, idx in enumerate(cluster_idxs)])

        if args.DBSCAN:
            avg_point = np.array([np.mean(cluster_pos[:, 0]), np.mean(cluster_pos[:, 1])])
            ins_seq = _pick_ins_seq(sub_dict)
            ins_len = _pick_ins_len(sub_dict)
            representing_var = make_representing_variant(
                variant, chrA, chrB,
                int(avg_point[0]), np.amin(cluster_pos[:, 0]), np.amax(cluster_pos[:, 0]),
                int(avg_point[1]), np.amin(cluster_pos[:, 1]), np.amax(cluster_pos[:, 1]),
                ins_seq, ins_len)
            cluster = [representing_var, sub_dict]
            f.write(vcf_line(cluster, f"cluster_{i}", sample_IDs, strip_chr) + "\n")
            i += 1
        else:
            i = overlap_cluster(sub_dict, sub_coords, variant, chrA, chrB, sample_IDs, args, f, i)

    return i


def export(args, sample_IDs):
    args.ins_seq_similarity = ins_similarity.resolve_ins_seq_threshold(args)
    with database.DB(args.db, memory=args.memory) as db:
        groups = db.query('SELECT DISTINCT var, chrA, chrB FROM SVDB')

        if any("INS" in var for var, _, _ in groups) and not db.has_ins_table():
            logger.warning(
                "database does not contain insertion sequence/length data — "
                "exporting insertions without sequence in ALT column. "
                "To enable full insertion export, run: svdb --build --upgrade --files <original_vcfs>"
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

    with open(args.prefix + ".vcf", 'w') as f:
        f.write(db_header(args) + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(sample_IDs)))
    export(args, sample_IDs)
