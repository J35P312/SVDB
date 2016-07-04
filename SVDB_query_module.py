import sys, os, glob
import readVCF
import SVDB_overlap_module
from operator import itemgetter
import SVDB_merge_vcf_module_cython

def main(args):
    #start by loading the variations
    variations = args.query_vcf
    queries = []
    if args.prefix:
        f=open(args.prefix+"_query.vcf","w")
    with open(variations) as fin:
        noOCCTag=1;
        infoFound=0;
        outputSource="Not_specified"
        for line in fin:
            #process and output the metadata
            if line.startswith("#") or line.startswith("="):
                #find the output source(cnvnator or Findtranslocations)
                meta_line=line.replace("#","");
                content=meta_line.split("=");
                if(content[0] == "source"):
                    outputSource=content[1].rstrip().split()[0];

                lookForFilter=meta_line.split("=");
                #the last infotag will be the Feature tag
                if(lookForFilter[0] != "INFO" and noOCCTag and infoFound==1):
                    if not args.prefix:
                        sys.stdout.write("##INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n".format(args.hit_tag));
                        sys.stdout.write("##INFO=<ID={},Number=1,Type=Float,Description=\"The frequency of the event in the database\">\n".format(args.frequency_tag));
                        sys.stdout.write(line);
                    else:
                        f.write("##INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n".format(args.hit_tag))
                        f.write("##INFO=<ID={},Number=1,Type=Float,Description=\"The frequency of the event in the database\">\n".format(args.frequency_tag))
                        f.write(line)

                    infoFound=0;noFeatureTag=0;
                elif(lookForFilter[0] == "INFO"):
                    if not args.prefix:
                        sys.stdout.write(line);
                    else:
                        f.write(line)

                    infoFound=1;
                    #there should only be one feature tag per vf file
                    if(line == "INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">".format(args.hit_tag)):
                        noOCCTag=0
                else:
                    if not args.prefix:
                        sys.stdout.write(line)
                    else:
                        f.write(line)
            else:
                #in this case I need to store a query
                chrA,posA,chrB,posB,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
                current_variation = [chrA, int(posA), chrB, int(posB),event_type, FORMAT, line] # plus a counter and the variation
                queries.append(current_variation)
                
    # at this point queries contains an entry for each variation
    #now query each sample.db present in the given folder and store the occurences
    db_file=args.db
    DBvariants={}
    db_size=1
    with open(db_file) as DB:
        for line in DB:
            if(line[0] != "#"):
                chrA,posA,chrB,posB,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
                if not chrA in DBvariants:
                    DBvariants[chrA]={}
                if not chrB in DBvariants[chrA]:
                    DBvariants[chrA][chrB]=[]
                info_field=line.split("\t")[7]
                DBvariants[chrA][chrB].append([int(posA),int(posB),event_type,FORMAT])
                db_size=len(FORMAT["GT"])

    
    for query in queries:
        hits = isVariationInDB(DBvariants, query,args)
        query[5] = hits
    for query in sorted(queries, key=itemgetter(5),reverse=args.invert):
        vcf_entry = query[6].rstrip()
        content=vcf_entry.split("\t")
        content[7]="{};{}={};{}={}".format(content[7],args.hit_tag, query[5],args.frequency_tag,(query[5]/float(db_size ) ))
        if not args.prefix:
            print(("\t").join(content))
        else:
            f.write(("\t").join(content)+"\n")


def isVariationInDB(DBvariants, Query_variant,args):
    chrA = Query_variant[0]
    chrApos = Query_variant[1]
    chrB =Query_variant[2]
    chrBpos = Query_variant[3]
    variation_type=Query_variant[4]
    samples=set([])
    if chrA in DBvariants:
        # now look if chrB is here
        if chrB in DBvariants[chrA]:
            # check if this variation is already present
            for event in DBvariants[chrA][chrB]:
                #check if the variant type of the events is the same
                if(event[2] == variation_type or args.no_var):
                    hit_tmp = None
                    if not args.ci:
                        if not (chrA == chrB):
                            hit_tmp=SVDB_overlap_module.precise_overlap(chrApos,chrBpos,event[0],event[1],args.bnd_distance)

                        elif chrBpos >= event[0] and event[1] >= chrApos:
                            hit_tmp = SVDB_overlap_module.isSameVariation(chrApos,chrBpos,event[0],event[1],args.overlap,args.bnd_distance)
                    else:
                        ciA_query,ciB_query,ciA_db,ciB_db=SVDB_merge_vcf_module_cython.find_ci(Query_variant,Query_variant)
                        hit_tmp=SVDB_overlap_module.ci_overlap(chrApos,chrBpos,ciA_query,ciB_query,event[0],event[1],[0,0],[0,0])

                    if hit_tmp != None:
                        genotype=event[-1]["GT"]
                        for i in range(0,len(genotype)):
                            GT=genotype[i]
                            if not GT == "0|0" and not GT == "0/0": 
                                samples = samples | set([i])
    hits=len(samples)
    return hits
