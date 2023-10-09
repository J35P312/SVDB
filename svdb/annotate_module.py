from __future__ import absolute_import

import gzip
import sys

import numpy as np

from . import database, overlap_module, readVCF


def main(args, output_file=None):
    # start by loading the variations
    queries = []
    if args.prefix:
        f = open(output_file, "w")
    infoFound = 0

    # at this point queries contains an entry for each variation
    # now query each sample.db present in the given folder and store the occurences
    db_file = args.db
    DBvariants = {}

    if not args.out_tags:
        args.out_tags=args.in_tags
    in_tag_to_out_tag_dictionary={}

    for i in range(len(args.in_tags)):
       in_tag_to_out_tag_dictionary[args.in_tags[i]]=args.out_tags[i]

    if not len(args.out_tags) == len(args.in_tags):
        print("Error: mismatching numbers of tags")
        print(args.out_tags,args.in_tags)
        quit()

    opener = gzip.open if db_file.endswith(".gz") else open
    # print FRQ_tag

    header_lines=[]
    in_tag_set=set(args.in_tags)

    with opener(db_file, 'rt') as lines:
       for line in lines:
           if line.startswith('#'):
              if "##INFO=<ID=" in line:
                 tag=line.split("##INFO=<ID=")[-1].split(",")[0] 
                 if tag in in_tag_set:
                    header_lines.append(  line.replace("##INFO=<ID={},".format(tag) ,"##INFO=<ID={},".format(in_tag_to_out_tag_dictionary[tag]))  )
              continue

           chrA, posA, chrB, posB, event_type, INFO, FORMAT = readVCF.readVCFLine(line)

           if chrA not in DBvariants:
               DBvariants[chrA] = {}
           if chrB not in DBvariants[chrA]:
               DBvariants[chrA][chrB] = {}
           if event_type not in DBvariants[chrA][chrB]:
               DBvariants[chrA][chrB][event_type] = {}
               DBvariants[chrA][chrB][event_type]["annotation"] = []
               DBvariants[chrA][chrB][event_type]["coordinates"] = []

           DBvariants[chrA][chrB][event_type]["coordinates"].append(np.array([int(posA), int(posB)]))

           annotation=[]
           for tag in range(len(args.in_tags)):
               if not args.in_tags[tag] in INFO:
                  continue

               annotation.append( "{}={}".format( args.out_tags[tag] ,INFO[args.in_tags[tag]] ) )

           DBvariants[chrA][chrB][event_type]["annotation"].append(annotation)


    for chrA in DBvariants:
        for chrB in DBvariants[chrA]:
            for var in DBvariants[chrA][chrB]:
                DBvariants[chrA][chrB][var]["coordinates"] = np.array(
                    DBvariants[chrA][chrB][var]["coordinates"])

    opener = gzip.open if args.query_vcf.endswith(".gz") else open
    writer = f.write if args.prefix else sys.stdout.write

    with opener(args.query_vcf, "rt") as lines:
        for line in lines:
            if line.startswith("#"):
                meta_line = line.replace("#", "")
                content = meta_line.split("=")

                lookForFilter = meta_line.split("=")
                # the last infotag will be the Feature tag
                if lookForFilter[0] != "INFO" and infoFound == 1:
                    for h in header_lines:
                        writer(h)

                    writer(line)
                    infoFound = 0
                elif lookForFilter[0] == "INFO":
                    writer(line)
                    infoFound = 1
                else:
                    if line[1] != "#":
                        writer("".join(header_lines))
                        writer("##SVDB_version={} cmd=\"{}\"\n".format(args.version, " ".join(sys.argv)))
                        writer(line)
                    else:
                        writer(line)
                continue

            # in this case I need to store a query
            chrA, posA, chrB, posB, event_type, INFO, FORMAT = readVCF.readVCFLine(line)
            # plus a counter and the variation
            queries.append([chrA, int(posA), chrB, int(posB), event_type, FORMAT, line])

    for query in queries:
        hits = queryVCFDB(DBvariants, query, args)
        query[5] = hits

    for query in queries:
        vcf_entry = query[6].strip()
        content = vcf_entry.split("\t")
        if query[5]:
            content[7] = ";".join([content[7],query[5][0] ])

        if not args.prefix:
            print(("\t").join(content))
        else:
            f.write(("\t").join(content) + "\n")

def queryVCFDB(DBvariants, query_variant, args):
    chrA = query_variant[0]
    chrApos = query_variant[1]
    chrB = query_variant[2]
    chrBpos = query_variant[3]
    variation_type = query_variant[4]

    annotations=[]
    similarity = []
    if chrA not in DBvariants:
            return 0
    if chrB not in DBvariants[chrA]:
            return 0

    for var in DBvariants[chrA][chrB]:
        if not args.no_var and variation_type != var:
            continue

        candidates = np.where((args.bnd_distance >= abs(DBvariants[chrA][chrB][var]["coordinates"][:, 0] - chrApos)) & (
            args.bnd_distance >= abs(DBvariants[chrA][chrB][var]["coordinates"][:, 1] - chrBpos)))
        if not len(candidates[0]) and not args.no_var:
                return 0

        # check if this variation is already present
        for candidate in candidates[0]:
            event = DBvariants[chrA][chrB][var]["coordinates"][candidate]
            sample_list = DBvariants[chrA][chrB][var]["annotation"][candidate]

            # check if the variant type of the events is the same
            hit_tmp = None
            match = False

            if not (chrA == chrB):
                hit_tmp, match = overlap_module.precise_overlap(
                    chrApos, chrBpos, event[0], event[1], args.bnd_distance)
            elif "INS" in variation_type:
                #insertions are treated as single points, overlap is not defined, and the maximum distance is determined by ins_distance
                hit_tmp, match = overlap_module.precise_overlap(
                    chrApos, chrBpos, event[0], event[1], args.ins_distance)
            else:
                hit_tmp, match = overlap_module.isSameVariation(
                    chrApos, chrBpos, event[0], event[1], args.overlap, args.bnd_distance)

            if match and sample_list:
                similarity.append(hit_tmp)
                annotations.append(sample_list)

    if match and annotations:
        if not (chrA == chrB):
           idx = similarity.index(min(similarity))
        else:
           idx = similarity.index(max(similarity))
        hits = annotations[idx]
    else:
        hits = 0

    return hits
