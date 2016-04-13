import sys, os, glob
import readVCF
import SVDB_overlap_module
from operator import itemgetter

def main(args):
    #start by loading the variations
    variations = args.query_vcf
    queries = []
    
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
                    sys.stdout.write("##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n");
                    sys.stdout.write("##INFO=<ID=FRQ,Number=1,Type=Float,Description=\"The frequency of the event in the database\">\n");
                    sys.stdout.write(line);
                    infoFound=0;noFeatureTag=0;
                elif(lookForFilter[0] == "INFO"):
                    sys.stdout.write(line);
                    infoFound=1;
                    #there should only be one feature tag per vf file
                    if(line == "INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">"):
                        noOCCTag=0
                else:
                    sys.stdout.write(line)
            else:
                #in this case I need to store a query
                chrA,posA,chrB,posB,event_type =readVCF.readVCFLine(line);
                current_variation = [chrA, int(posA), chrB, int(posB),event_type, 0, line] # plus a counter and the variation
                queries.append(current_variation)
                
    # at this point queries contains an entry for each variation
    #now query each sample.db present in the given folder and store the occurences
    db_file=args.db
    DBvariants={}
    with open(db_file) as DB:
        for line in DB:
            if(line[0] != "#"):
                chrA,posA,chrB,posB,event_type =readVCF.readVCFLine(line);
                if not chrA in DBvariants:
                    DBvariants[chrA]={}
                if not chrB in DBvariants[chrA]:
                    DBvariants[chrA][chrB]=[]
                info_field=line.split("\t")[7]
                DBvariants[chrA][chrB].append([int(posA),int(posB),event_type,info_field])

    db_size=1
    for query in queries:
        hits,db = isVariationInDB(DBvariants, query,args)
        query[5] = hits
        if db != 0:
            db_size=db

    for query in sorted(queries, key=itemgetter(5)):
        vcf_entry = query[6].rstrip()
        content=vcf_entry.split("\t")
        content[7]="{};OCC={};FRQ={}".format(content[7], query[5],(query[5]/float(db_size ) ))
        print(("\t").join(content))



def isVariationInDB(DBvariants, Query_variant,args):
    chrA = Query_variant[0]
    chrApos = Query_variant[1]
    chrB =Query_variant[2]
    chrBpos = Query_variant[3]
    variation_type=Query_variant[4]

    db_size=0;
    samples=set([])
    if chrA in DBvariants:
        # now look if chrB is here
        if chrB in DBvariants[chrA]:
            # check if this variation is already present
            for event in DBvariants[chrA][chrB]:
                #check if the variant type of the events is the same
                if(event[2] == variation_type or args.no_var):
                
                    if not (chrA == chrB):
                        hit_tmp=SVDB_overlap_module.precise_overlap(chrApos,chrBpos,event[0],event[1],args.bnd_distance)
                        if hit_tmp != None:
                            hit_tag=event[-1].strip().split(";SAMPLES=")[-1]
                            db_tag=event[-1].split("NSAMPLES=")[-1]
                            
                            hit_tag=hit_tag.split(";")[0];
                            db_tag=db_tag.split(";")[0];
                            
                            for hit in hit_tag.split("|"):
                                samples.add(hit)
                            db_size=int(db_tag)
                    elif chrBpos >= event[0] and event[1] >= chrApos:
                        hit_tmp = SVDB_overlap_module.isSameVariation(chrApos,chrBpos,event[0],event[1],args.overlap)
                        if hit_tmp != None:
                            hit_tag=event[-1].split(";SAMPLES=")[-1]
                            db_tag=event[-1].split("NSAMPLES=")[-1]
                            
                            hit_tag=hit_tag.split(";")[0];
                            db_tag=db_tag.split(";")[0];
                            
                            for hit in hit_tag.split("|"):
                                samples.add(hit)
                            db_size=int(db_tag)
    hits=len(samples);
    return hits,db_size
