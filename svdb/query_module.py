from __future__ import absolute_import
import sys, os, glob
from . import readVCF
from . import overlap_module
from operator import itemgetter
from . import merge_vcf_module_cython
import sqlite3

import numpy as np

def main(args):
    #start by loading the variations
    queries = []
    if args.prefix:
        f=open(args.prefix+"_query.vcf","w")
    noOCCTag=1;
    infoFound=0;
        
    for line in open(args.query_vcf):
        if line[0] == "#":
            meta_line=line.replace("#","");
            content=meta_line.split("=");

            lookForFilter=meta_line.split("=");
            #the last infotag will be the Feature tag
            if(lookForFilter[0] != "INFO" and noOCCTag and infoFound==1):
                if not args.prefix:
                    sys.stdout.write("##INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n".format(args.out_occ));
                    sys.stdout.write("##INFO=<ID={},Number=1,Type=Float,Description=\"The frequency of the event in the database\">\n".format(args.out_frq));
                    sys.stdout.write(line);
                else:
                    f.write("##INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n".format(args.out_occ))
                    f.write("##INFO=<ID={},Number=1,Type=Float,Description=\"The frequency of the event in the database\">\n".format(args.out_frq))
                    f.write(line)

                infoFound=0;noFeatureTag=0;
            elif(lookForFilter[0] == "INFO"):
                if not args.prefix:
                    sys.stdout.write(line);
                else:
                    f.write(line)

                infoFound=1;
                #there should only be one feature tag per vf file
                if(line == "INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">".format(args.out_occ)):
                    noOCCTag=0


            else:
                if line[1] != "#":
                    if not args.prefix:
                       print "##SVDB_version={} cmd=\"{}\"".format(args.version," ".join(sys.argv))
                    else:
                       f.write( "##SVDB_version={} cmd=\"{}\"".format(args.version," ".join(sys.argv)) )
                if not args.prefix:
                    sys.stdout.write(line)
                else:
                    f.write(line)
            continue

        #in this case I need to store a query
        chrA,posA,chrB,posB,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
        current_variation = [chrA, int(posA), chrB, int(posB),event_type, FORMAT, line] # plus a counter and the variation
        queries.append(current_variation)
                
    # at this point queries contains an entry for each variation
    #now query each sample.db present in the given folder and store the occurences
    
    if args.bedpedb or args.db:
        if args.bedpedb:
            args.db=args.bedpedb
        db_file=args.db
        DBvariants={}
        db_size=1
        Use_OCC_tag=False
        if args.in_occ:
            OCC_tag=args.in_occ
            Use_OCC_tag=True           

        if args.in_frq:
            FRQ_tag=args.in_frq

        #print FRQ_tag
        for line in open(db_file):
            if line[0] == "#":
                continue
            
            if args.bedpedb:
                content=line.strip().split()

                if ( content[0] == content[2] and (int(content[1]) < int(content[3])) ) or (content[0] < content[2]):
                     chrA=content[0]
                     posA=int(content[1])
                     chrB=content[2]
                     posB=int(content[3])
                else:
                     chrA=content[2]
                     posA=int(content[3])
                     chrB=content[0]
                     posB=int(content[1])

                event_type=content[4]
                hits=int(content[5])
                frequency=float(content[6])
                FORMAT=[False]

            else:
                chrA,posA,chrB,posB,event_type,INFO,FORMAT = readVCF.readVCFLine(line)

            if not chrA in DBvariants:
                DBvariants[chrA]={}
            if not chrB in DBvariants[chrA]:
                DBvariants[chrA][chrB]={}
            if not event_type in DBvariants[chrA][chrB]:
                DBvariants[chrA][chrB][event_type]={}
                DBvariants[chrA][chrB][event_type]["samples"]=[]
                DBvariants[chrA][chrB][event_type]["coordinates"]=[]

            DBvariants[chrA][chrB][event_type]["coordinates"].append(np.array([int(posA),int(posB)]))
            if "GT" in FORMAT and not Use_OCC_tag:
                DBvariants[chrA][chrB][event_type]["samples"].append(np.array(FORMAT["GT"]))
                db_size=len(FORMAT["GT"])
            elif args.bedpedb:
                DBvariants[chrA][chrB][event_type]["samples"].append([hits,frequency])
                Use_OCC_tag=True
   
            else:
                try:
                    OCC=INFO[OCC_tag]
                    FRQ=INFO[FRQ_tag]
                    DBvariants[chrA][chrB][event_type]["samples"].append([OCC,FRQ])
                    Use_OCC_tag=True
                except:
                    print("Error: frequency or hit tag not found! Make sure to set the --in_occ AND --in_frq to the number and frequency of alleles/individuals as presented in the INFO column of the input db\n" )
                    print("database variants not having the --in_occ or --in_frq tag must be removed")
                    print("you may also skip these parameters and cluster based on the GT entry of the format column (if such exists)" )
                    quit()

        for chrA in DBvariants:
            for chrB in DBvariants[chrA]:
                for var in DBvariants[chrA][chrB]:
                    DBvariants[chrA][chrB][var]["coordinates"]=np.array(DBvariants[chrA][chrB][var]["coordinates"])                
                    DBvariants[chrA][chrB][var]["samples"]=np.array(DBvariants[chrA][chrB][var]["samples"])
                 
        for query in queries:
            hits = queryVCFDB(DBvariants, query,args,Use_OCC_tag)
            query[5] = hits

        for query in queries:
            vcf_entry = query[6].strip()
            content=vcf_entry.split("\t")
            if not Use_OCC_tag:
                if query[5]:
                    content[7]="{};{}={};{}={}".format(content[7],args.out_occ, query[5],args.out_frq,(query[5]/float(db_size ) ))
            else:
                if query[5][0]:
                    content[7]="{};{}={};{}={}".format(content[7],args.out_occ, int(query[5][0]),args.out_frq,query[5][1])    
            
            if not args.prefix:
                print(("\t").join(content))
            else:
                f.write(("\t").join(content)+"\n")
        return()
            
                
    elif args.sqdb:
        db_file=args.sqdb
        conn = sqlite3.connect(args.sqdb)
        
        if args.memory:
            memory_db=sqlite3.connect(':memory:')
            db_dump="".join(line for line in conn.iterdump())
            memory_db.executescript(db_dump)
            conn.close()
            c = memory_db.cursor()
        else:
            c=conn.cursor()
    
        db_size=0
        A='SELECT DISTINCT sample FROM SVDB'
        for sample in c.execute(A):
            db_size +=1
        if not db_size:
            print "error: no samples in the db"
            quit()

        for query in queries:
            hits = SQDB(query,args,c)
            query[5] = hits    

            
            
            
    for query in queries:
        vcf_entry = query[6].strip()
        content=vcf_entry.split("\t")
        if query[5]:
            content[7]="{};{}={};{}={}".format(content[7],args.out_occ, query[5],args.out_frq,(query[5]/float(db_size ) ))
        if not args.prefix:
            print(("\t").join(content))
        else:
            f.write(("\t").join(content)+"\n")



def queryVCFDB(DBvariants, Query_variant,args,Use_OCC_tag):
    chrA = Query_variant[0]
    chrApos = Query_variant[1]
    chrB =Query_variant[2]
    chrBpos = Query_variant[3]
    variation_type=Query_variant[4]
    samples=set([])
    frequency=[]
    occ=[]
    similarity=[]
    if not chrA in DBvariants:
        if Use_OCC_tag:
            return([0,0])
        else:
            return 0
    if not chrB in DBvariants[chrA]:
        if Use_OCC_tag:
            return([0,0])
        else:
            return 0
    for var in DBvariants[chrA][chrB]:
        if not args.no_var and variation_type != var:
            continue

        #candidates=DBvariants[chrA][chrB][var]["coordinates"][ ( args.bnd_distance >= abs(DBvariants[chrA][chrB][var]["coordinates"][:,0] - chrApos)  ) & ( args.bnd_distance >= abs(DBvariants[chrA][chrB][var]["coordinates"][:,1] - chrBpos)  ) ]
        candidates=np.where( ( args.bnd_distance >= abs(DBvariants[chrA][chrB][var]["coordinates"][:,0] - chrApos)  ) & ( args.bnd_distance >= abs(DBvariants[chrA][chrB][var]["coordinates"][:,1] - chrBpos)  ) ) 
        if not len(candidates[0]) and not args.no_var:
            if Use_OCC_tag:
                return([0,0])
            else:
                return 0
        # check if this variation is already present
        for candidate in candidates[0]:
            event=DBvariants[chrA][chrB][var]["coordinates"][candidate]
            sample_list=DBvariants[chrA][chrB][var]["samples"][candidate]
            #check if the variant type of the events is the same
            hit_tmp = None
            match=False

            if not (chrA == chrB):
                hit_tmp,match=overlap_module.precise_overlap(chrApos,chrBpos,event[0],event[1],args.bnd_distance)
            elif chrBpos >= event[0] and event[1] >= chrApos:
                hit_tmp,match = overlap_module.isSameVariation(chrApos,chrBpos,event[0],event[1],args.overlap,args.bnd_distance)

            if match:
                similarity.append(hit_tmp)
                if Use_OCC_tag:
                    occ.append(sample_list[0])
                    frequency.append(sample_list[1])
                else:
                    for i in range(0,len(sample_list)):
                        GT=sample_list[i]
                        if not GT == "0|0" and not GT == "0/0": 
                            samples = samples | set([i])
    if Use_OCC_tag:
        if occ:
            if not (chrA == chrB):
                idx=similarity.index(min(similarity))
            else:
                idx=similarity.index(max(similarity))
            hits=[ occ[idx],frequency[idx] ]
        else:
            hits=[0,0]
    else:
        hits=len(samples)
    
    return hits
    

def SQDB(Query_variant,args,c):
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
    occurances=0
    for hit in hits:
        if not ci and variant["chrA"] == variant["chrB"]:
            var={}
            var["posA"]=int( hit[0] )
            var["posB"]=int( hit[1] )
            var["index"]=hit[2]
            similar=overlap_module.isSameVariation(variant["posA"],variant["posB"],var["posA"],var["posB"],overlap,distance)
            if similar:
                match.add( var["index"] )
                    
        else:
            match.add(hit[0])

    occurances=len(match)
    return occurances  
