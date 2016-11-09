import sys, os, glob
import readVCF
import SVDB_overlap_module
from operator import itemgetter
import SVDB_merge_vcf_module_cython
import sqlite3

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
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
            continue

        #in this case I need to store a query
        chrA,posA,chrB,posB,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
        current_variation = [chrA, int(posA), chrB, int(posB),event_type, FORMAT, line] # plus a counter and the variation
        queries.append(current_variation)
                
    # at this point queries contains an entry for each variation
    #now query each sample.db present in the given folder and store the occurences
    
    if args.db:
        db_file=args.db
        DBvariants={}
        db_size=1
        for line in open(db_file):
            if line[0] == "#":
                continue
            
            chrA,posA,chrB,posB,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
            if not chrA in DBvariants:
                DBvariants[chrA]={}
            if not chrB in DBvariants[chrA]:
                DBvariants[chrA][chrB]=[]
            info_field=line.split("\t")[7]
            DBvariants[chrA][chrB].append([int(posA),int(posB),event_type,FORMAT])
            db_size=len(FORMAT["GT"])
        for query in queries:
            hits = queryVCFDB(DBvariants, query,args)
            query[5] = hits            
                
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

        if not args.DBSCAN:
            for query in queries:
                hits = SQDB(query,args,c)
                query[5] = hits    
        else:
            queries= DBSCAN_query(queries,args,c)
            
            
            
    for query in sorted(queries, key=itemgetter(5),reverse=args.invert):
        vcf_entry = query[6].rstrip()
        content=vcf_entry.split("\t")
        content[7]="{};{}={};{}={}".format(content[7],args.hit_tag, query[5],args.frequency_tag,(query[5]/float(db_size ) ))
        if not args.prefix:
            print(("\t").join(content))
        else:
            f.write(("\t").join(content)+"\n")



def queryVCFDB(DBvariants, Query_variant,args):
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
            similar=SVDB_overlap_module.isSameVariation(variant["posA"],variant["posB"],var["posA"],var["posB"],overlap,distance)
            if similar:
                match.add( var["index"] )
                    
        else:
            match.add(hit[0])

    occurances=len(match)
    return occurances
    
def DBSCAN_query(queries,args,c):
    query_db={}
    for i in range(0,len(queries)):
        query = queries[i]
        
        if not query[0] in query_db:
            query_db[query[0]] ={}
        if not query[2] in query_db[query[0]]:
            query_db[query[0]][query[2]] ={}
        if not query[4] in query_db[query[0]][query[2]]:
            query_db[query[0]][query[2]][query[4]] = []          
        
        query_db[query[0]][query[2]][query[4]].append([query[1],query[3], i,0])

    for chrA in query_db:
        for chrB in query_db[chrA]:
            for variant in query_db[chrA][chrB]:
                query_db[chrA][chrB][variant] = np.array(query_db[chrA][chrB][variant])

    chrA_list=[]    
    for chrA in c.execute('SELECT DISTINCT chrA FROM SVDB'):
        chrA_list.append(chrA[0])
        
    chrB_list=[]    
    for chrB in c.execute('SELECT DISTINCT chrB FROM SVDB'):
        chrB_list.append(chrB[0])
    
    for chrA in chrA_list:
        if not chrA in query_db:
            continue
        for chrB in chrB_list:
            if not chrB in query_db[chrA]:
                continue
                   
            chr_db={}            
            for hit in c.execute('SELECT posA,posB,sample,idx,var FROM SVDB WHERE chrA == \'{}\' AND chrB == \'{}\''.format(chrA,chrB)):
                if not hit[-1] in chr_db:
                    chr_db[ hit[-1] ]={}
                    chr_db[ hit[-1] ]["coordinates"]=[]
                    chr_db[ hit[-1] ]["var_info"]=[]
                    chr_db[ hit[-1] ]["index"]=[]
                
                chr_db[ hit[-1] ]["coordinates"].append([hit[0],hit[1],hit[2],-1])
        
            for variant in chr_db:
                if not variant in query_db[chrA][chrB]:
                    continue            
                db_size=len(chr_db[variant]["coordinates"])
                chr_db[variant]["coordinates"]=np.array(chr_db[variant]["coordinates"])
                chr_db[variant]["var_info"]=np.array(chr_db[variant]["var_info"])
                chr_db[variant]["index"]=np.array(chr_db[variant]["index"])
                  
                variants=np.vstack( (chr_db[variant]["coordinates"],query_db[chrA][chrB][variant]) )                  
                db = DBSCAN(eps=args.epsilon, min_samples=args.min_pts).fit(variants[:,0:2])
                
                core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                core_samples_mask[db.core_sample_indices_] = True
                labels = db.labels_

                stats=[]
                query_labels=labels[db_size:]
                query_variants=variants[db_size:]
                for i in range(0,len(query_db[chrA][chrB][variant])):
                    variant_cluster= query_labels[i]
                    if variant_cluster == -1:
                        continue
                    else:
                        #compute the number of samples having this particular variant
                        db_label_list=labels[0:db_size]
                        query_db[chrA][chrB][variant][i][-1] = len( set(   variants[ np.where(db_label_list == variant_cluster) ][:,2] ) )
                        if( len( set(   variants[ np.where(db_label_list == variant_cluster) ][:,2] ) ) == 0 ):
                            print "ERROR"

    for chrA in query_db:
        for chrB in query_db[chrA]:
            for variant in query_db[chrA][chrB]:
                for entry in query_db[chrA][chrB][variant]:
                    queries[ entry[-2] ][5] = entry[-1]
    return(queries)        
