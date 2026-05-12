import logging
import sys

import numpy as np

from . import database, ins_similarity, overlap_module, read_vcf, vcf_utils

logger = logging.getLogger(__name__)

# Sequence gate applies only within this positional hard cap, regardless of ins_distance.
_INS_SEQ_HARD_CAP = 25


def _read_query_vcf(args, writer):
    """Read the query VCF, rewrite the header to writer, return collected variant queries."""
    queries = []
    noOCCTag = 1
    infoFound = 0
    db_path = args.db or args.bedpedb or args.sqdb

    with vcf_utils.open_vcf(args.query_vcf) as lines:
        for line in lines:
            if line.startswith("#"):
                meta_line = line.replace("#", "")
                lookForFilter = meta_line.split("=")
                # inject OCC/FRQ INFO tags after the last INFO line
                if lookForFilter[0] != "INFO" and noOCCTag and infoFound == 1:
                    writer(
                        f'##INFO=<ID={args.out_occ},Number=1,Type=Integer,Description="The number of occurrences of the event in the database {db_path}">\n'
                    )
                    writer(
                        f'##INFO=<ID={args.out_frq},Number=1,Type=Float,Description="The frequency of the event in the database {db_path}">\n'
                    )
                    writer(line)
                    infoFound = 0
                elif lookForFilter[0] == "INFO":
                    writer(line)
                    infoFound = 1
                    if (
                        line
                        == f'INFO=<ID={args.out_occ},Number=1,Type=Integer,Description="The number of occurrences of the event in the database {db_path}">'
                    ):
                        noOCCTag = 0
                else:
                    if line[1] != "#":
                        writer(
                            f'##SVDB_version={args.version} cmd="{" ".join(sys.argv)}"\n'
                        )
                    writer(line)
                continue

            variant = read_vcf.readVCFLine(line)
            if variant is None:
                logger.debug("skipping unparseable line: %s", line.rstrip())
                continue
            queries.append([variant.chrA, int(variant.posA), variant.chrB, int(variant.posB), variant.event_type, variant.fmt, line, variant.ins_seq, variant.svlen])

    return queries


def _load_vcf_db(args):
    """Load a VCF or BEDPE database into an in-memory lookup structure.

    Returns (DBvariants, db_size, use_OCC_tag).
    """
    if args.bedpedb:
        args.db = args.bedpedb
    db_file = args.db
    DBvariants = {}
    db_size = 1
    use_OCC_tag = False
    OCC_tag = None
    FRQ_tag = None

    if args.in_occ:
        OCC_tag = args.in_occ
        use_OCC_tag = True
    if args.in_frq:
        FRQ_tag = args.in_frq

    with vcf_utils.open_vcf(db_file) as lines:
        for line in lines:
            if line.startswith('#'):
                continue

            if args.bedpedb:
                content = line.strip().split()
                if (content[0] == content[2] and (int(content[1]) < int(content[3]))) or (content[0] < content[2]):
                    chrA = content[0]
                    posA = int(content[1])
                    chrB = content[2]
                    posB = int(content[3])
                else:
                    chrA = content[2]
                    posA = int(content[3])
                    chrB = content[0]
                    posB = int(content[1])
                event_type = content[4]
                hits = int(content[5])
                frequency = float(content[6])
                FORMAT = [False]
                INFO = {}
            else:
                v = read_vcf.readVCFLine(line)
                chrA, posA, chrB, posB, event_type, INFO, FORMAT = v.chrA, v.posA, v.chrB, v.posB, v.event_type, v.info, v.fmt

            if chrA not in DBvariants:
                DBvariants[chrA] = {}
            if chrB not in DBvariants[chrA]:
                DBvariants[chrA][chrB] = {}
            if event_type not in DBvariants[chrA][chrB]:
                DBvariants[chrA][chrB][event_type] = {
                    "samples": [], "coordinates": [], "svlens": [], "sequences": []}

            DBvariants[chrA][chrB][event_type]["coordinates"].append(np.array([int(posA), int(posB)]))

            # Store SVLEN and sequence for insertion entries (VCF format only, not BEDPE)
            if not args.bedpedb and "INS" in event_type:
                svlen = v.svlen
                seq = v.ins_seq
            else:
                svlen = None
                seq = ""
            DBvariants[chrA][chrB][event_type]["svlens"].append(svlen)
            DBvariants[chrA][chrB][event_type]["sequences"].append(seq)
            if "GT" in FORMAT and not use_OCC_tag:
                DBvariants[chrA][chrB][event_type]["samples"].append(np.array(FORMAT["GT"]))
                db_size = len(FORMAT["GT"])
            elif args.bedpedb:
                DBvariants[chrA][chrB][event_type]["samples"].append([hits, frequency])
                use_OCC_tag = True
            else:
                try:
                    OCC = INFO[OCC_tag]
                    FRQ = INFO[FRQ_tag]
                    DBvariants[chrA][chrB][event_type]["samples"].append([OCC, FRQ])
                    use_OCC_tag = True
                except KeyError:
                    logger.warning(
                        "skipping db variant at %s:%s — OCC/FRQ tag not found in INFO "
                        "(expected %s / %s); use --in_occ/--in_frq to set the correct tag names",
                        chrA, posA, OCC_tag, FRQ_tag,
                    )
                    # remove entries added above for this variant to keep all lists in sync
                    DBvariants[chrA][chrB][event_type]["coordinates"].pop()
                    DBvariants[chrA][chrB][event_type]["svlens"].pop()
                    DBvariants[chrA][chrB][event_type]["sequences"].pop()
                    continue

    for chrA in DBvariants:
        for chrB in DBvariants[chrA]:
            for var in DBvariants[chrA][chrB]:
                DBvariants[chrA][chrB][var]["coordinates"] = np.array(DBvariants[chrA][chrB][var]["coordinates"])
                DBvariants[chrA][chrB][var]["samples"] = np.array(DBvariants[chrA][chrB][var]["samples"])

    return DBvariants, db_size, use_OCC_tag


def _write_vcfdb_results(queries, args, writer, db_size, use_OCC_tag):
    """Write annotated query results for VCF/BEDPE database queries."""
    for query in queries:
        vcf_entry = query[6].strip()
        content = vcf_entry.split("\t")
        if not use_OCC_tag:
            if query[5]:
                content[7] = f"{content[7]};{args.out_occ}={query[5]};{args.out_frq}={query[5] / float(db_size)}"
        else:
            if query[5][0]:
                content[7] = f"{content[7]};{args.out_occ}={int(query[5][0])};{args.out_frq}={query[5][1]}"
        writer(("\t").join(content) + "\n")


def _write_sqdb_results(queries, args, writer, db_size):
    """Write annotated query results for SQLite database queries."""
    for query in queries:
        vcf_entry = query[6].strip()
        content = vcf_entry.split("\t")
        frq = query[5] / float(db_size)
        if frq > args.max_frq:
            continue
        if query[5]:
            content[7] = f"{content[7]};{args.out_occ}={query[5]};{args.out_frq}={frq}"
        writer(("\t").join(content) + "\n")


def main(args, output_file=None):
    # Resolve insertion sequence similarity threshold before any matching.
    args.ins_seq_similarity = ins_similarity.resolve_ins_seq_threshold(args)

    if args.prefix:
        f = open(output_file, "w")
        writer = f.write
    else:
        writer = sys.stdout.write

    queries = _read_query_vcf(args, writer)

    if args.bedpedb or args.db:
        DBvariants, db_size, use_OCC_tag = _load_vcf_db(args)

        for query in queries:
            query[5] = queryVCFDB(DBvariants, query, args, use_OCC_tag)

        _write_vcfdb_results(queries, args, writer, db_size, use_OCC_tag)
        return None

    elif args.sqdb:
        db = database.DB(db=args.sqdb, memory=args.memory)
        db_size = len(db)
        if not db_size:
            logger.error("no samples found in the database")
            sys.exit(1)

        has_ins_queries = any("INS" in q[4] for q in queries)
        has_ins_table = db.has_ins_table()
        if has_ins_queries and not has_ins_table:
            logger.warning(
                "database does not contain insertion sequence/length data — "
                "matching on position only. To enable full insertion matching, "
                "run: svdb --build --upgrade --files <original_vcfs>"
            )

        for query in queries:
            query[5] = SQDB(query, args, db, has_ins_table)

        _write_sqdb_results(queries, args, writer, db_size)


def queryVCFDB(DBvariants, query_variant, args, use_OCC_tag):
    chrA = query_variant[0]
    chrApos = query_variant[1]
    chrB = query_variant[2]
    chrBpos = query_variant[3]
    variation_type = query_variant[4]
    samples = set([])
    frequency = []
    occ = []
    similarity = []

    is_ins = "INS" in variation_type
    query_seq = query_variant[7] if is_ins else ""
    query_svlen = query_variant[8] if is_ins else None

    ins_svlen_ratio = getattr(args, "ins_svlen_ratio", 0.90)
    ins_seq_threshold = getattr(args, "ins_seq_similarity", 0.75)
    no_ins_seq = getattr(args, "no_ins_seq", False)
    if chrA not in DBvariants:
        if use_OCC_tag:
            return([0, 0])
        else:
            return 0
    if chrB not in DBvariants[chrA]:
        if use_OCC_tag:
            return([0, 0])
        else:
            return 0
    for var in DBvariants[chrA][chrB]:
        if not args.no_var and variation_type != var:
            continue

        candidates = np.where((args.bnd_distance >= abs(DBvariants[chrA][chrB][var]["coordinates"][:, 0] - chrApos)) & (
            args.bnd_distance >= abs(DBvariants[chrA][chrB][var]["coordinates"][:, 1] - chrBpos)))
        if not len(candidates[0]) and not args.no_var:
            if use_OCC_tag:
                return([0, 0])
            else:
                return 0
        # check if this variation is already present
        for candidate in candidates[0]:
            event = DBvariants[chrA][chrB][var]["coordinates"][candidate]
            sample_list = DBvariants[chrA][chrB][var]["samples"][candidate]
            hit_tmp = None
            match = False

            if not (chrA == chrB):
                hit_tmp, match = overlap_module.precise_overlap(
                    chrApos, chrBpos, event[0], event[1], args.bnd_distance)
            elif is_ins:
                # insertions are treated as single points; max distance determined by ins_distance
                hit_tmp, match = overlap_module.precise_overlap(
                    chrApos, chrBpos, event[0], event[1], args.ins_distance)
            else:
                hit_tmp, match = overlap_module.isSameVariation(
                    chrApos, chrBpos, event[0], event[1], args.overlap, args.bnd_distance)

            # SVLEN ratio gate (VCF db mode only — sequences/svlens stored at load time)
            if match and is_ins and chrA == chrB:
                db_svlen = DBvariants[chrA][chrB][var]["svlens"][candidate]
                if db_svlen is not None and query_svlen is not None:
                    if not overlap_module.insertion_svlen_match(query_svlen, db_svlen, ins_svlen_ratio):
                        match = False

            # Sequence similarity gate (within hard cap; skip if no_ins_seq or no sequences)
            if match and is_ins and chrA == chrB and not no_ins_seq:
                pos_dist = abs(chrApos - int(event[0]))
                if pos_dist <= _INS_SEQ_HARD_CAP:
                    db_seq = DBvariants[chrA][chrB][var]["sequences"][candidate]
                    if not ins_similarity.sequence_gate(query_seq, db_seq, ins_seq_threshold):
                        match = False

            if match:
                similarity.append(hit_tmp)
                if use_OCC_tag:
                    occ.append(sample_list[0])
                    frequency.append(sample_list[1])
                else:
                    for i in range(0, len(sample_list)):
                        GT = sample_list[i]
                        if not GT == "0|0" and not GT == "0/0":
                            samples = samples | set([i])
    if use_OCC_tag:
        if occ:
            if not (chrA == chrB):
                idx = similarity.index(min(similarity))
            elif "INS" in variation_type:
                idx = similarity.index(min(similarity))
            else:
                idx = similarity.index(max(similarity))
            hits = [occ[idx], frequency[idx]]
        else:
            hits = [0, 0]
    else:
        hits = len(samples)

    return hits


def SQDB(query_variant, args, db, has_ins_table=False):
    is_ins = "INS" in query_variant[4]
    distance = getattr(args, "ins_distance", 25) if is_ins else args.bnd_distance
    overlap = args.overlap
    variant = {"type": query_variant[4],
               "chrA": query_variant[0], "posA": query_variant[1],
               "chrB": query_variant[2], "posB": query_variant[3]}

    use_ins_table = is_ins and has_ins_table and not getattr(args, "no_ins_seq", False)

    if use_ins_table:
        selection = "s.posA, s.posB, s.sample, s.idx, i.ins_seq, i.ins_len"
        join = "LEFT JOIN INS i ON s.idx = i.idx"
        table = "SVDB s"
        A = (
            f"SELECT {selection} FROM {table} {join} "
            f"WHERE s.var == '{variant['type']}' AND s.chrA == '{variant['chrA']}' "
            f"AND s.chrB == '{variant['chrB']}' "
            f"AND s.posA <= {variant['posA'] + distance} AND s.posA >= {variant['posA'] - distance} "
            f"AND s.posB <= {variant['posB'] + distance} AND s.posB >= {variant['posB'] - distance}"
        )
    else:
        selection = "posA, posB, sample" if variant["chrA"] == variant["chrB"] else "sample"
        A = (
            f"SELECT {selection} FROM SVDB "
            f"WHERE var == '{variant['type']}' AND chrA == '{variant['chrA']}' "
            f"AND chrB == '{variant['chrB']}' "
            f"AND posA <= {variant['posA'] + distance} AND posA >= {variant['posA'] - distance} "
            f"AND posB <= {variant['posB'] + distance} AND posB >= {variant['posB'] - distance}"
        )

    hits = db.query(A)

    ins_svlen_ratio = getattr(args, "ins_svlen_ratio", 0.90)
    ins_seq_threshold = getattr(args, "ins_seq_similarity", 0.75)
    query_seq = query_variant[7] if is_ins else ""
    query_svlen = query_variant[8] if is_ins else None

    match = set()
    for hit in hits:
        if variant["chrA"] == variant["chrB"]:
            hit_posA, hit_posB, hit_sample = int(hit[0]), int(hit[1]), hit[2]
            if use_ins_table:
                hit_idx, hit_seq, hit_len = hit[3], hit[4], hit[5]
                _, similar = overlap_module.precise_overlap(
                    variant["posA"], variant["posB"], hit_posA, hit_posB, distance)
                if similar and hit_len is not None and query_svlen is not None:
                    if not overlap_module.insertion_svlen_match(query_svlen, hit_len, ins_svlen_ratio):
                        similar = False
                pos_dist = abs(variant["posA"] - hit_posA)
                if similar and pos_dist <= _INS_SEQ_HARD_CAP:
                    if not ins_similarity.sequence_gate(query_seq, hit_seq or "", ins_seq_threshold):
                        similar = False
                if similar:
                    match.add(hit_idx)
            elif is_ins:
                _, similar = overlap_module.precise_overlap(
                    variant["posA"], variant["posB"], hit_posA, hit_posB, distance)
                if similar:
                    match.add(hit_sample)
            else:
                similar, _ = overlap_module.isSameVariation(
                    variant["posA"], variant["posB"], hit_posA, hit_posB, overlap, distance)
                if similar:
                    match.add(hit_sample)
        else:
            match.add(hit[0])

    return len(match)
