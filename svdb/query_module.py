import logging
import sys

import numpy as np

from . import database, overlap_module, read_vcf, vcf_utils

logger = logging.getLogger(__name__)


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
            queries.append([variant.chrA, int(variant.posA), variant.chrB, int(variant.posB), variant.event_type, variant.fmt, line])

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
                DBvariants[chrA][chrB][event_type] = {"samples": [], "coordinates": []}

            DBvariants[chrA][chrB][event_type]["coordinates"].append(np.array([int(posA), int(posB)]))
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
                    logger.error(
                        "frequency or hit tag not found — set --in_occ and --in_frq to the "
                        "OCC/FRQ tags in the INFO column of the input db (usually AC/OCC and AF/FRQ)"
                    )
                    logger.error("database variants without those tags must be removed first")
                    logger.error(
                        "alternatively omit --in_occ/--in_frq and let SVDB use the GT field instead"
                    )
                    sys.exit(1)

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
                content[7] = "{};{}={};{}={}".format(content[7], args.out_occ, query[5], args.out_frq, (query[5] / float(db_size)))
        else:
            if query[5][0]:
                content[7] = "{};{}={};{}={}".format(content[7], args.out_occ, int(query[5][0]), args.out_frq, query[5][1])
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

        for query in queries:
            query[5] = SQDB(query, args, db)

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
            elif "INS" in variation_type:
                # insertions are treated as single points; max distance determined by ins_distance
                hit_tmp, match = overlap_module.precise_overlap(
                    chrApos, chrBpos, event[0], event[1], args.ins_distance)
            else:
                hit_tmp, match = overlap_module.isSameVariation(
                    chrApos, chrBpos, event[0], event[1], args.overlap, args.bnd_distance)

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


def SQDB(query_variant, args, db):
    distance = args.bnd_distance
    overlap = args.overlap
    variant = {"type": query_variant[4],
               "chrA": query_variant[0], "posA": query_variant[1], "ci_A_start": 0, "ci_A_end": 0,
               "chrB": query_variant[2], "posB": query_variant[3], "ci_B_start": 0, "ci_B_end": 0}

    selection = "posA, posB, sample" if variant["chrA"] == variant["chrB"] else "sample"

    A = 'SELECT {} FROM SVDB WHERE var == \'{}\' AND chrA == \'{}\' AND chrB == \'{}\' AND posA <= {} AND posA >= {} AND posB <= {} AND posB >= {}'.format(
        selection, variant["type"], variant["chrA"], variant["chrB"],
        variant["posA"] + distance, variant["posA"] - distance,
        variant["posB"] + distance, variant["posB"] - distance)
    hits = db.query(A)

    match = set([])
    for hit in hits:
        if variant["chrA"] == variant["chrB"]:
            var = {"posA": int(hit[0]), "posB": int(hit[1]), "index": hit[2]}
            similar, _ = overlap_module.isSameVariation(variant["posA"], variant["posB"],
                                                        var["posA"], var["posB"], overlap, distance)
            if similar:
                match.add(var["index"])
        else:
            match.add(hit[0])

    return len(match)
