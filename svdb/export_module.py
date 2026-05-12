import logging
import sys
from collections import Counter

import numpy as np

from . import dbscan, database, ins_similarity, overlap_module

logger = logging.getLogger(__name__)


def make_representing_variant(variant_type, chrA, chrB, posA, ci_A_start, ci_A_end, posB, ci_B_start, ci_B_end, ins_seq=None):
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
    }


def build_genotype_columns(sample_IDs, hit_sample_ids):
    """Return per-sample GT strings in sample_IDs order.

    Samples present in hit_sample_ids get './1'; all others get '0/0'.
    """
    zygosity = {sample: "0/0" for sample in sample_IDs}
    for sid in hit_sample_ids:
        zygosity[sid] = "./1"
    return [zygosity[sample] for sample in sample_IDs]


def fetch_index_variant(db, index):
    has_ins = db.has_ins_table()
    if has_ins:
        A = (
            'SELECT s.posA, s.ci_A_lower, s.ci_A_upper, s.posB, s.ci_B_lower, s.ci_B_upper, '
            's.sample, i.ins_seq, i.ins_len '
            'FROM SVDB s LEFT JOIN INS i ON s.idx = i.idx '
            'WHERE s.idx IN ({})'.format(", ".join([str(idx) for idx in index]))
        )
    else:
        A = 'SELECT posA, ci_A_lower, ci_A_upper, posB, ci_B_lower, ci_B_upper, sample FROM SVDB WHERE idx IN ({}) '.format(
            ", ".join([str(idx) for idx in index]))
    hits = db.query(A)
    variant = {}
    coordinates = []
    for i, hit in enumerate(hits):
        variant[i] = {
            "posA": int(hit[0]),
            "ci_A_start": int(hit[1]),
            "ci_A_end": int(hit[2]),
            "posB": int(hit[3]),
            "ci_B_start": int(hit[4]),
            "ci_B_end": int(hit[5]),
            "sample_id": hit[6],
            "ins_seq": hit[7] if has_ins else None,
            "ins_len": hit[8] if has_ins else None,
        }
        coordinates.append([i, int(hit[0]), int(hit[3])])
    return variant, np.array(coordinates)


def fetch_cluster_variant(db, index):
    has_ins = db.has_ins_table()
    if has_ins:
        query = (
            'SELECT s.posA, s.posB, s.sample, s.idx, i.ins_seq, i.ins_len '
            'FROM SVDB s LEFT JOIN INS i ON s.idx = i.idx '
            'WHERE s.idx IN ({})'.format(", ".join([str(idx) for idx in index]))
        )
    else:
        query = 'SELECT posA, posB, sample, idx FROM SVDB WHERE idx IN ({}) '.format(
                ", ".join([str(idx) for idx in index]))
    hits = db.query(query)

    variant_dict = {}
    for hit in hits:
        i = int(hit[3])
        variant_dict[i] = {
            "posA": int(hit[0]),
            "posB": int(hit[1]),
            "sample_id": hit[2],
            "ins_seq": hit[4] if has_ins else None,
            "ins_len": hit[5] if has_ins else None,
        }
    return variant_dict


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


def vcf_line(cluster, id_tag, sample_IDs):
    info_field = "SVTYPE={};".format(cluster[0]["type"])
    vcf_line = []
    vcf_line.append(cluster[0]["chrA"])
    vcf_line.append(str(cluster[0]["posA"]))
    vcf_line.append(id_tag)
    vcf_line.append("N")
    is_ins = cluster[0]["type"] == "INS" and cluster[0]["chrA"] == cluster[0]["chrB"]
    if cluster[0]["chrA"] == cluster[0]["chrB"] and cluster[0]["type"] != "BND":
        ins_seq = cluster[0].get("ins_seq")
        if is_ins and ins_seq:
            vcf_line.append("N" + ins_seq)
            info_field += "SVLEN={};".format(len(ins_seq))
        else:
            vcf_line.append("<" + cluster[0]["type"] + ">")
            info_field += "END={};SVLEN={};".format(cluster[0]["posB"], abs(cluster[0]["posA"] - cluster[0]["posB"]))
    else:
        vcf_line.append("N[{}:{}[".format(cluster[0]["chrB"], cluster[0]["posB"]))

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


def fetch_variants(variant, chrA, chrB, db):
    chr_db = {}
    chr_db[variant] = {}

    hits = db.query(f'SELECT posA,posB,sample,idx,var FROM SVDB WHERE var == \'{variant}\'AND chrA == \'{chrA}\' AND chrB == \'{chrB}\'')
    if not hits:
        return False

    x = [v[0] for v in hits]
    y = [v[1] for v in hits]

    chr_db[variant]["coordinates"] = np.column_stack((x, y))
    chr_db[variant]["var_info"] = np.array([v[2] for v in hits])
    chr_db[variant]["index"] = np.array([v[3] for v in hits])
    return chr_db


def _pick_ins_seq(variant_dict):
    """Return the most common non-null ins_seq across the cluster, or None."""
    seqs = [v.get("ins_seq") for v in variant_dict.values() if v.get("ins_seq")]
    if not seqs:
        return None
    return Counter(seqs).most_common(1)[0][0]


def overlap_cluster(db, indexes, variant, chrA, chrB, sample_IDs, args, f, i):
    variant_dictionary, coordinates = fetch_index_variant(db, indexes)
    if "INS" in variant:
        similarity_matrix = expand_chain(
            variant_dictionary, coordinates, chrA, chrB, args.ins_distance, -1,
            ins_svlen_ratio=getattr(args, "ins_svlen_ratio", None),
            ins_seq_threshold=getattr(args, "ins_seq_similarity", None),
            no_ins_seq=getattr(args, "no_ins_seq", False))
    else:
        similarity_matrix = expand_chain(
           variant_dictionary, coordinates, chrA, chrB, args.bnd_distance, args.overlap)

    clusters = cluster_variants(variant_dictionary, similarity_matrix)
    for clustered_variants in clusters:
        clustered_variants[0]["type"] = variant
        clustered_variants[0]["chrA"] = chrA
        clustered_variants[0]["chrB"] = chrB
        clustered_variants[0]["ins_seq"] = _pick_ins_seq(clustered_variants[1])
        f.write(vcf_line(clustered_variants, f"cluster_{i}", sample_IDs) + "\n")
    return i + len(clusters)


def svdb_cluster_main(chrA, chrB, variant, sample_IDs, args, db, i, f):
    chr_db = fetch_variants(variant, chrA, chrB, db)
    if not chr_db:
        return i

    #DBSCAN clustering according to the user set parameters
    if args.DBSCAN:
        labels = dbscan.main(chr_db[variant]["coordinates"], args.epsilon, args.min_pts)
    elif "INS" in variant:
        #insertions are clustered based on the ins_distance, which is typically smaller than the BND_distance
        labels = dbscan.main(chr_db[variant]["coordinates"], args.ins_distance, 2)
    else:
        #clustering of all other variants
        labels = dbscan.main(chr_db[variant]["coordinates"], args.bnd_distance, 2)

    unique_labels = set(labels)
    # print the unique variants
    unique_xy = chr_db[variant]["coordinates"][labels == -1]
    unique_index = chr_db[variant]["index"][labels == -1]
    for xy, indexes in zip(unique_xy, unique_index):
        variant_dictionary = fetch_cluster_variant(db, [indexes])
        ins_seq = _pick_ins_seq(variant_dictionary)
        representing_var = make_representing_variant(
            variant, chrA, chrB, xy[0], xy[0], xy[0], xy[1], xy[1], xy[1], ins_seq)
        cluster = [representing_var, variant_dictionary]
        f.write(vcf_line(cluster, f"cluster_{i}", sample_IDs) + "\n")
        i += 1
    del unique_xy
    del unique_index

    # print the clusters
    for unique_label in unique_labels:
        if unique_label == -1:
            continue
        class_member_mask = (labels == unique_label)
        xy = chr_db[variant]["coordinates"][class_member_mask]
        indexes = chr_db[variant]["index"][class_member_mask]

        if args.DBSCAN:
            avg_point = np.array([np.mean(xy[:, 0]), np.mean(xy[:, 1])])

            variant_dictionary = fetch_cluster_variant(db, indexes)
            ins_seq = _pick_ins_seq(variant_dictionary)

            representing_var = make_representing_variant(
                variant, chrA, chrB,
                int(avg_point[0]), np.amin(xy[:, 0]), np.amax(xy[:, 0]),
                int(avg_point[1]), np.amin(xy[:, 1]), np.amax(xy[:, 1]),
                ins_seq)
            cluster = [representing_var, variant_dictionary]
            f.write(vcf_line(cluster, f"cluster_{i}", sample_IDs) + "\n")
            i += 1

        else:
            i = overlap_cluster(db, indexes, variant, chrA,
                                chrB, sample_IDs, args, f, i)

    return i


def export(args, sample_IDs):
    args.ins_seq_similarity = ins_similarity.resolve_ins_seq_threshold(args)
    with database.DB(args.db, memory=args.memory) as db:
        chrA_list = db.query_column('SELECT DISTINCT chrA FROM SVDB')
        chrB_list = db.query_column('SELECT DISTINCT chrB FROM SVDB')
        var_list = db.query_column('SELECT DISTINCT var FROM SVDB')

        if any("INS" in v for v in var_list) and not db.has_ins_table():
            logger.warning(
                "database does not contain insertion sequence/length data — "
                "exporting insertions without sequence in ALT column. "
                "To enable full insertion export, run: svdb --build --upgrade --files <original_vcfs>"
            )

        i = 0
        with open(args.prefix + ".vcf", 'a') as f:
            for chrA in chrA_list:
                for chrB in chrB_list:
                    for variant in var_list:
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
