from __future__ import absolute_import

import gzip
import sys

import numpy as np

from db import DB

from . import overlap_module, readVCF


def main(args):
    # start by loading the variations
    queries = []
    if args.prefix:
        f = open(args.prefix + "_query.vcf", "w")
    noOCCTag = 1
    infoFound = 0

    opener = gzip.open if args.query_vcf.endswith(".gz") else open
    writer = f.write if args.prefix else sys.stdout.write

    with opener(args.query_vcf, "rt") as lines:
        for line in lines:
            if line.startswith("#"):
                meta_line = line.replace("#", "")
                content = meta_line.split("=")

                lookForFilter = meta_line.split("=")
                # the last infotag will be the Feature tag
                if lookForFilter[0] != "INFO" and noOCCTag and infoFound == 1:
                    writer("##INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n".format(args.out_occ))
                    writer("##INFO=<ID={},Number=1,Type=Float,Description=\"The frequency of the event in the database\">\n".format(args.out_frq))
                    writer(line)
                    infoFound = 0
                elif lookForFilter[0] == "INFO":
                    writer(line)
                    infoFound = 1
                    # there should only be one feature tag per vf file
                    if line == "INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">".format(args.out_occ):
                        noOCCTag = 0
                else:
                    if line[1] != "#":
                        writer("##SVDB_version={} cmd=\"{}\"".format(args.version, " ".join(sys.argv)))
                        writer(line)
                continue

            # in this case I need to store a query
            chrA, posA, chrB, posB, event_type, INFO, FORMAT = readVCF.readVCFLine(line)
            # plus a counter and the variation
            queries.append([chrA, int(posA), chrB, int(posB), event_type, FORMAT, line])

    # at this point queries contains an entry for each variation
    # now query each sample.db present in the given folder and store the occurences

    if args.bedpedb or args.db:
        if args.bedpedb:
            args.db = args.bedpedb
        db_file = args.db
        DBvariants = {}
        db_size = 1
        use_OCC_tag = False
        if args.in_occ:
            OCC_tag = args.in_occ
            use_OCC_tag = True

        if args.in_frq:
            FRQ_tag = args.in_frq

        opener = gzip.open if db_file.endswith(".gz") else open
        # print FRQ_tag
        with opener(db_file, 'rt') as lines:
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

                else:
                    chrA, posA, chrB, posB, event_type, INFO, FORMAT = readVCF.readVCFLine(line)

                if chrA not in DBvariants:
                    DBvariants[chrA] = {}
                if chrB not in DBvariants[chrA]:
                    DBvariants[chrA][chrB] = {}
                if event_type not in DBvariants[chrA][chrB]:
                    DBvariants[chrA][chrB][event_type] = {}
                    DBvariants[chrA][chrB][event_type]["samples"] = []
                    DBvariants[chrA][chrB][event_type]["coordinates"] = []

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
                        print("Error: frequency or hit tag not found! Make sure to set the --in_occ AND --in_frq to the number and frequency of alleles/individuals as presented in the INFO column of the input db\n")
                        print("database variants not having the --in_occ or --in_frq tag must be removed")
                        print("you may also skip these parameters and cluster based on the GT entry of the format column (if such exists)")
                        quit()

        for chrA in DBvariants:
            for chrB in DBvariants[chrA]:
                for var in DBvariants[chrA][chrB]:
                    DBvariants[chrA][chrB][var]["coordinates"] = np.array(
                        DBvariants[chrA][chrB][var]["coordinates"])
                    DBvariants[chrA][chrB][var]["samples"] = np.array(
                        DBvariants[chrA][chrB][var]["samples"])

        for query in queries:
            hits = queryVCFDB(DBvariants, query, args, use_OCC_tag)
            query[5] = hits

        for query in queries:
            vcf_entry = query[6].strip()
            content = vcf_entry.split("\t")
            if not use_OCC_tag:
                if query[5]:
                    content[7] = "{};{}={};{}={}".format(content[7], args.out_occ, query[5], args.out_frq, (query[5] / float(db_size)))
            else:
                if query[5][0]:
                    content[7] = "{};{}={};{}={}".format(content[7], args.out_occ, int(query[5][0]), args.out_frq, query[5][1])

            if not args.prefix:
                print(("\t").join(content))
            else:
                f.write(("\t").join(content) + "\n")
        return None

    elif args.sqdb:
        db = DB(db=args.sqdb, memory=args.memory)

        db_size = len(db)
        if not db_size:
            # TODO: Raise a custom DB exception
            print("error: no samples in the db")
            quit()

        for query in queries:
            hits = SQDB(query, args, db)
            query[5] = hits

    for query in queries:
        vcf_entry = query[6].strip()
        content = vcf_entry.split("\t")
        frq = query[5] / float(db_size)
        if frq > args.max_frq:
            continue
        if query[5]:
            content[7] = "{};{}={};{}={}".format(content[7], args.out_occ, query[5], args.out_frq, frq)
        if not args.prefix:
            print(("\t").join(content))
        else:
            f.write(("\t").join(content) + "\n")


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
            # check if the variant type of the events is the same
            hit_tmp = None
            match = False

            if not (chrA == chrB):
                hit_tmp, match = overlap_module.precise_overlap(
                    chrApos, chrBpos, event[0], event[1], args.bnd_distance)
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
