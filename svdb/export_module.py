from __future__ import absolute_import
import argparse
from . import readVCF
from . import overlap_module
import subprocess
import sqlite3
import glob
import os
import sys
from . import DBSCAN

import numpy as np
import time


def fetch_index_variant(c,index):
    A='SELECT posA ,ci_A_lower ,ci_A_upper ,posB ,ci_B_lower ,ci_B_upper, sample FROM SVDB WHERE idx IN ({}) '.format( ", ".join([ str(idx) for idx in index ]) )
    hits = c.execute(A) 
    variant={}
    coordinates=[]
    i=0
    for hit in hits:
        variant[i]={}
        variant[i]["posA"]=int( hit[0] )
        variant[i]["ci_A_start"]= int(hit[1])
        variant[i]["ci_A_end"]= int(hit[2])
        variant[i]["posB"]=int( hit[3] )
        variant[i]["ci_B_start"]= int(hit[4])
        variant[i]["ci_B_end"]= int(hit[5])
        variant[i]["sample_id"]= hit[6]
        coordinates.append([i,int(hit[0]),int(hit[3])])
        i += 1

    return(variant, np.array(coordinates) )

def fetch_cluster_variant(c,index):
    A='SELECT posA, posB, sample, idx FROM SVDB WHERE idx IN ({}) '.format( ", ".join([ str(idx) for idx in index ]) )
    hits = c.execute(A)
    
    variant_dict={}
    for hit in hits:
        variant_dict[ int( hit[3] ) ]={}
        variant_dict[ int( hit[3] ) ]["posA"]=int( hit[0] )
        variant_dict[ int( hit[3] ) ]["posB"]=int( hit[1] )
        variant_dict[ int( hit[3] ) ]["sample_id"]=hit[2]
    return(variant_dict)


def db_header(args):
    headerString="##fileformat=VCFv4.1\n";
    headerString+="##source=SVDB\n";
    headerString+="##ALT=<ID=DEL,Description=\"Deletion\">\n";
    headerString+="##ALT=<ID=DUP,Description=\"Duplication\">\n";
    headerString+="##ALT=<ID=INV,Description=\"Inversion\">\n";
    headerString+="##ALT=<ID=BND,Description=\"Break end\">\n";
    headerString+="##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    headerString+="##INFO=<ID=END,Number=1,Type=String,Description=\"End of an intra-chromosomal variant\">\n";
    headerString+="##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurences of the event in the database\">\n";
    headerString+="##INFO=<ID=NSAMPLES,Number=1,Type=Integer,Description=\"the number of samples within the database\">\n";
    headerString+="##INFO=<ID=VARIANTS,Number=1,Type=Integer,Description=\"a| separated list of the positions of the clustered variants\">\n";
    headerString+="##INFO=<ID=FRQ,Number=1,Type=Float,Description=\"the frequency of the variant\">\n";
    headerString+="##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"
    headerString+="##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n"
    headerString+="##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">\n"
    headerString+="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    headerString+="##SVDB_version={} cmd=\"{}\"".format(args.version," ".join(sys.argv))
    return(headerString)


def vcf_line(cluster,id_tag,sample_IDs):
    info_field="SVTYPE={};".format(cluster[0]["type"])
    vcf_line=[]
    vcf_line.append( cluster[0]["chrA"] )
    vcf_line.append( str(cluster[0]["posA"]) )
    vcf_line.append( id_tag )
    vcf_line.append( "N" )
    if cluster[0]["chrA"] == cluster[0]["chrB"] and not cluster[0]["type"] == "BND":
        vcf_line.append( "<" + cluster[0]["type"] + ">"  )
        info_field += "END={};SVLEN={};".format(cluster[0]["posB"],abs(cluster[0]["posA"]-cluster[0]["posB"]))
    else:
        vcf_line.append( "N[{}:{}[".format(cluster[0]["chrB"],cluster[0]["posB"]) )
        
        
    sample_set=set([])
    CIPOS=[]
    CIEND=[]

    for variant in cluster[1]:
        CIPOS.append(cluster[1][variant]["posA"])
        CIEND.append(cluster[1][variant]["posB"])
        sample_set.add(cluster[1][variant]["sample_id"])

    CIPOS_start=-abs(cluster[0]["posA"]-min(CIPOS))
    CIPOS_end=abs(cluster[0]["posA"]-max(CIPOS))

    CIEND_start=-abs(cluster[0]["posB"]-min(CIEND))
    CIEND_end=abs(cluster[0]["posB"]-max(CIEND))

    info_field += "NSAMPLES={};OCC={};FRQ={};CIPOS={},{};CIEND={},{};".format(len(sample_IDs),len(sample_set),round(len(sample_set)/float(len(sample_IDs)) ,4),CIPOS_start,CIPOS_end,CIEND_start,CIEND_end)
    variant_field="VARIANTS="
    for variant in cluster[1]:
        variant_field +="|{}:{}:{}".format(cluster[1][variant]["sample_id"],cluster[1][variant]["posA"],cluster[1][variant]["posB"])
    info_field+=variant_field
    vcf_line.append( "." )
    vcf_line.append( "PASS" )
    vcf_line.append( info_field )
    zygosity_list={}
    for sample in sample_IDs:
        zygosity_list[sample]="0/0"
    
    for variant in cluster[1]:
        zygosity_list[cluster[1][variant]["sample_id"]]="./1"
    format=[]
    for sample in sample_IDs:
        format.append(zygosity_list[sample])
    vcf_line.append("GT")
    vcf_line.append("\t".join(format))
    return( "\t".join(vcf_line) )

def expand_chain(chain,coordinates,chrA,chrB,distance,overlap):
    chain_data={}
    for i in range(0,len(chain)):
        variant=chain[i]
        chain_data[i]=[]

        if chrA == chrB:
            rows=coordinates[ ( distance >= abs(coordinates[:,1] - variant["posA"])  ) & ( distance >= abs(coordinates[:,2] - variant["posB"])  ) & ( variant["posB"] >=  coordinates[:,1] )  & (coordinates[:,2] >= variant["posB"] ) ]

        else:
            rows=coordinates[ ( distance >= abs(coordinates[:,1] - variant["posA"])  ) & ( distance >= abs(coordinates[:,2] - variant["posB"])  ) ]
        candidates= rows[:,0]
        for j in candidates:

            var=chain[j]
            similar = False
            if chrA != chrB:
                similar=True
            else:
                similar,match=overlap_module.isSameVariation(variant["posA"],variant["posB"],var["posA"],var["posB"],overlap,distance)
            if match:
                chain_data[i].append(j)
        chain_data[i]=np.array(chain_data[i])
    return(chain_data)

def cluster_variants(variant_dictionary,similarity_matrix):
    
    cluster_sizes=[]
    for i in range(0,len(variant_dictionary)):
        cluster_sizes.append([i,len(similarity_matrix[i])])
    
    clusters=[]
    for index in sorted(cluster_sizes,key=lambda x: (x[1]), reverse=True):
        i=index[0]
        if similarity_matrix[i][0] == -1:
            continue

        cluster_dictionary={}
        for var in similarity_matrix[i]:
            similarity_matrix[var][0] = -1
            cluster_dictionary[var]=variant_dictionary[var]
        variant=variant_dictionary[i]
       
        clusters.append( [variant,cluster_dictionary] )      
    return(clusters)

def fetch_variants(variant,chrA,chrB,c):
    chr_db={}
    chr_db[ variant ]={}
           
    hits = c.execute('SELECT posA,posB,sample,idx,var FROM SVDB WHERE var == \'{}\'AND chrA == \'{}\' AND chrB == \'{}\''.format(variant,chrA,chrB)).fetchall()
    if not hits:
        return False

    x=[v[0] for v in hits]
    y=[v[1] for v in hits]

    chr_db[variant]["coordinates"]=np.column_stack((x,y))
    chr_db[variant]["var_info"]=np.array([v[2] for v in hits])
    chr_db[variant]["index"]=np.array([v[3] for v in hits])
    return chr_db


def overlap_cluster(c,indexes,variant,chrA,chrB,sample_IDs,args,f,i):
    variant_dictionary,coordinates=fetch_index_variant(c,indexes)
    similarity_matrix=expand_chain(variant_dictionary,coordinates,chrA,chrB,args.bnd_distance,args.overlap)
    clusters=cluster_variants(variant_dictionary,similarity_matrix)
    for clustered_variants in clusters:
        clustered_variants[0]["type"]=variant
        clustered_variants[0]["chrA"]=chrA
        clustered_variants[0]["chrB"]=chrB
        f.write( vcf_line(clustered_variants,"cluster_{}".format(i),sample_IDs )+"\n")
        i += 1
    return i

def svdb_cluster_main(chrA,chrB,variant,sample_IDs,args,c,i):
    f = open(args.prefix+".vcf",'a')

    chr_db=fetch_variants(variant,chrA,chrB,c)
    if not chr_db:
        f.close()
        return i

    if args.DBSCAN:
        db=DBSCAN.main(chr_db[variant]["coordinates"],args.epsilon,args.min_pts)
    else:
        db=DBSCAN.main(chr_db[variant]["coordinates"],args.bnd_distance,2)

    unique_labels = set(db)
    stats=[]
    #print the unique variants
    unique_xy=chr_db[variant]["coordinates"][  db == -1 ]
    unique_index=chr_db[variant]["index"][ db == -1 ]
    for j in range(0,len(unique_xy) ):
        xy = unique_xy[j]
        indexes=[ unique_index[j] ]
                    
                    
        variant_dictionary=fetch_cluster_variant(c,indexes)
                        
        representing_var={}
        representing_var["type"]=variant
        representing_var["chrA"]=chrA
        representing_var["chrB"]=chrB
        representing_var["posA"]= xy[0]
        representing_var["ci_A_start"]=  xy[0]
        representing_var["ci_A_end"]= xy[0]
        representing_var["posB"]= xy[1]
        representing_var["ci_B_start"]= xy[1]
        representing_var["ci_B_end"]= xy[1]                      
                        

        cluster=[representing_var,variant_dictionary] 
        f.write( vcf_line(cluster,"cluster_{}".format(i),sample_IDs )+"\n")
        i += 1
    del unique_xy
    del unique_index
                      
    #print the clusters
    for k in unique_labels:
        if k == -1:
            continue
        class_member_mask = (db == k)
        xy = chr_db[variant]["coordinates"][class_member_mask]
        indexes=chr_db[variant]["index"][class_member_mask]

        if args.DBSCAN:
             avg_point=np.array([np.mean(xy[:,0]),np.mean(xy[:,1])])
                        
             distance=[]
             for point in xy:
                distance.append( ( sum((avg_point-point)**2))**0.5  )
                        
             variant_dictionary=fetch_cluster_variant(c,indexes)
 
             representing_var={}
             representing_var["type"]=variant
             representing_var["chrA"]=chrA
             representing_var["chrB"]=chrB
             representing_var["posA"]= int(avg_point[0])
             representing_var["ci_A_start"]=  np.amin(xy[:,0])
             representing_var["ci_A_end"]= np.amax(xy[:,0])
             representing_var["posB"]= int(avg_point[1])
             representing_var["ci_B_start"]= np.amin(xy[:,1])
             representing_var["ci_B_end"]= np.amax(xy[:,1])                      
                            
             cluster=[representing_var,variant_dictionary]   
             f.write( vcf_line(cluster,"cluster_{}".format(i),sample_IDs )+"\n")
             i += 1
        else:
            i=overlap_cluster(c,indexes,variant,chrA,chrB,sample_IDs,args,f,i)
    f.close()
    return i

def export(args,sample_IDs):

    db=args.db
    conn = sqlite3.connect(db)
    if args.memory:
        memory_db=sqlite3.connect(':memory:')
        db_dump="".join(line for line in conn.iterdump())
        memory_db.executescript(db_dump)
        conn.close
        c = memory_db.cursor()
    else:
        c=conn.cursor()

        
    chrA_list=[]    
    for chrA in c.execute('SELECT DISTINCT chrA FROM SVDB'):
        chrA_list.append(chrA[0])
        
    chrB_list=[]    
    for chrB in c.execute('SELECT DISTINCT chrB FROM SVDB'):
        chrB_list.append(chrB[0])

    var_list=[]
    for variant in c.execute('SELECT DISTINCT var FROM SVDB'):
        var_list.append(variant[0])

    i=0;
    for chrA in chrA_list:
        for chrB in chrB_list:
            for variant in var_list:
                i=svdb_cluster_main(chrA,chrB,variant,sample_IDs,args,c,i)

 
def main(args):
    sample_IDs=[]
    if not args.prefix:
        args.prefix=args.db.replace(".db","")
    
    
    conn = sqlite3.connect(args.db)
    c=conn.cursor()

    A='SELECT DISTINCT sample FROM SVDB'
    for sample in c.execute(A):
            sample_IDs.append(sample[0])

    conn.close()

    f = open(args.prefix+".vcf",'w')
    f.write( db_header(args)+"\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(sample_IDs)))
    f.close()
    export(args,sample_IDs)
