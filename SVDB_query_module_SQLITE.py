import sys, os, glob
import readVCF
import SVDB_overlap_module
from operator import itemgetter
import SVDB_merge_vcf_module_cython
import sqlite3

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
    memory_db=sqlite3.connect(':memory:')
    conn = sqlite3.connect(args.db)
    db_dump="".join(line for line in conn.iterdump())
    memory_db.executescript(db_dump)
    conn.close()
    c = memory_db.cursor()
    #c=conn.cursor()
    
    db_size=0
    A='SELECT DISTINCT sample FROM SVDB'
    for sample in c.execute(A):
            db_size +=1
    if not db_size:
        print "error: no samples in the db"
        quit()

    for query in queries:
        hits = isVariationInDB(query,args,c)
        query[5] = hits
    for query in sorted(queries, key=itemgetter(5),reverse=args.invert):
        vcf_entry = query[6].rstrip()
        content=vcf_entry.split("\t")
        content[7]="{};{}={};{}={}".format(content[7],args.hit_tag, query[5],args.frequency_tag,(query[5]/float(db_size ) ))
        if not args.prefix:
            print(("\t").join(content))
        else:
            f.write(("\t").join(content)+"\n")


def isVariationInDB(Query_variant,args,c):
    ci = args.ci
    distance = args.bnd_distance
    overlap = args.overlap
    variant={}
    variant["type"]=Query_variant[4]
    variant["chrA"]=Query_variant[0]
    variant["chrB"]=Query_variant[2]
    variant["posA"]=Query_variant[1]
    variant["ci_A_start"]= 0
    variant["ci_A_end"]= 0
    variant["posB"]=Query_variant[3]
    variant["ci_B_start"]= 0
    variant["ci_B_end"]= 0
            
    selection ="sample"
    if variant["chrA"] == variant["chrB"]:
        selection = "posA, posB, sample"
    if not ci:
        A='SELECT {} FROM SVDB WHERE var == \'{}\' AND chrA == \'{}\' AND chrB == \'{}\' AND posA <= {} AND posA >= {} AND posB <= {} AND posB >= {}'.format(selection,variant["type"],variant["chrA"],variant["chrB"],variant["posA"]+distance, variant["posA"] -distance,variant["posB"] + distance, variant["posB"]-distance)
                #print A
    else:
        A='SELECT idx FROM SVDB WHERE var == \'{}\' AND chrA == \'{}\' AND chrB == \'{}\' AND posA <= {} AND posA >= {} AND posB <= {} AND posB >= {}'.format(variant["type"],variant["chrA"],variant["chrB"],variant["posA"]+variant["ci_A_end"], variant["posA"] -variant["ci_A_start"],variant["posB"] + variant["ci_B_end"], variant["posB"] - variant["ci_B_start"])
            
    hits = c.execute(A)
    match=set([])
    for hit in hits:
        if not ci and variant["chrA"] == variant["chrB"]:
            var={}
            var["posA"]=int( hit[0] )
            var["posB"]=int( hit[1] )
            var["index"]=hit[2]
            similar=SVDB_overlap_module.isSameVariation(variant["posA"],variant["posB"],var["posA"],var["posB"],overlap,distance)
            if similar:
                match.add( var["index"] )
                    
        else:
            match.add(hit[0])

        hits=len(match)
    return hits
