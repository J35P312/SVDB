from __future__ import absolute_import
import argparse
from . import readVCF
from . import overlap_module
import subprocess
import sqlite3
import glob
import os

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import numpy as np



def fetch_index_variant(c,index):
    A='SELECT var ,chrA , chrB ,posA ,ci_A_lower ,ci_A_upper ,posB ,ci_B_lower ,ci_B_upper, sample, idx  FROM SVDB WHERE idx IN ({}) '.format( ", ".join([ str(idx) for idx in index ]) )
    hits = c.execute(A) 
    variant={}
    for hit in hits:
        variant[ int(hit[10])]={}
        variant[ int(hit[10])]["type"]=hit[0]
        variant[ int(hit[10])]["chrA"]=hit[1]
        variant[ int(hit[10])]["chrB"]=hit[2]
        variant[ int(hit[10])]["posA"]=int( hit[3] )
        variant[ int(hit[10])]["ci_A_start"]= int(hit[4])
        variant[ int(hit[10])]["ci_A_end"]= int(hit[5])
        variant[ int(hit[10])]["posB"]=int( hit[6] )
        variant[ int(hit[10])]["ci_B_start"]= int(hit[7])
        variant[ int(hit[10])]["ci_B_end"]= int(hit[8])
        variant[ int(hit[10])]["sample_id"]= hit[9]
    return(variant)

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


def db_header():
    headerString="##fileformat=VCFv4.1\n";
    headerString+="##source=SVDB\n";
    headerString+="##ALT=<ID=DEL,Description=\"Deletion\">\n";
    headerString+="##ALT=<ID=DUP,Description=\"Duplication\">\n";
    headerString+="##ALT=<ID=INV,Description=\"Inversion\">\n";
    headerString+="##ALT=<ID=BND,Description=\"Break end\">\n";
    headerString+="##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    headerString+="##INFO=<ID=END,Number=1,Type=String,Description=\"End of an intra-chromosomal variant\">\n";
    headerString+="##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n";
    headerString+="##INFO=<ID=NSAMPLES,Number=1,Type=Integer,Description=\"the number of samples within the database\">\n";
    headerString+="##INFO=<ID=VARIANTS,Number=1,Type=Integer,Description=\"a| separated list of the positions of the clustered variants\">\n";
    headerString+="##INFO=<ID=FRQ,Number=1,Type=Float,Description=\"the frequency of the vriant\">\n";
    headerString+="##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"
    headerString+="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    return(headerString)


def vcf_line(cluster,id_tag,sample_IDs):
    info_field="SVTYPE={};".format(cluster[0]["type"])
    vcf_line=[]
    vcf_line.append( cluster[0]["chrA"] )
    vcf_line.append( str(cluster[0]["posA"]) )
    vcf_line.append( id_tag )
    vcf_line.append( "N" )
    if cluster[0]["chrA"] == cluster[0]["chrB"]:
        vcf_line.append( "<" + cluster[0]["type"] + ">"  )
        info_field += "END={};SVLEN={};".format(cluster[0]["posB"],abs(cluster[0]["posA"]-cluster[0]["posB"]))
    else:
        vcf_line.append( "N[{}:{}[".format(cluster[0]["chrB"],cluster[0]["posB"]) )
        
        
    sample_set=set([])
    for variant in cluster[1]:
        sample_set.add(cluster[1][variant]["sample_id"])
    info_field += "NSAMPLES={};OCC={};FRQ={};".format(len(sample_IDs),len(sample_set),round(len(sample_set)/float(len(sample_IDs)) ,2))
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
        zygosity_list[cluster[1][variant]["sample_id"]]="./."
    format=[]
    for sample in sample_IDs:
        format.append(zygosity_list[sample])
    vcf_line.append("GT")
    vcf_line.append("\t".join(format))
    return( "\t".join(vcf_line) )

def expand_chain(chain,distance,overlap,ci):

    chain_data={}
    for index in chain:
        variant=chain[index]
        chain_data[index]=[]
        for var_index in chain:
            var=chain[var_index]

            similar = False              
            if ci and var["posA"] >= variant["ci_A_start"]  and variant["ci_A_end"] >= var["posA"] and var["posB"] >= variant["ci_B_start"] and variant["ci_B_end"] >= var["posB"]:
                similar=True
            else:
                if variant["chrA"] == variant["chrB"]:
                     similar=overlap_module.isSameVariation(variant["posA"],variant["posB"],var["posA"],var["posB"],overlap,distance)
                else:
                     similar=overlap_module.precise_overlap(variant["posA"],variant["posB"],var["posA"],var["posB"],distance)
            if similar:
                chain_data[index].append(var_index)

    return(chain_data)

def cluster_variants(variant_dictionary,chain,chain_data):
    clusters=[]
    for index in sorted(chain_data, key=lambda index: len(chain_data[index]), reverse=True):
        if not chain:
            break
        if not index in chain:
            continue

        cluster_dictionary={}
        for var in chain_data[index]:
            if var in chain:
                chain.remove(var)
            cluster_dictionary[var]=variant_dictionary[var]
        variant=variant_dictionary[index]
       
        clusters.append( [variant,cluster_dictionary] )      
    return(clusters)

def export(args,sample_IDs):

    db=args.db
    f = open(args.prefix+".vcf",'a')
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
                chr_db={}
                chr_db[ variant ]={}
            
                hits = c.execute('SELECT posA,posB,sample,idx,var FROM SVDB WHERE var == \'{}\'AND chrA == \'{}\' AND chrB == \'{}\''.format(variant,chrA,chrB)).fetchall()
                if not hits:
                   continue

                x=[v[0] for v in hits]
                y=[v[1] for v in hits]

                chr_db[variant]["coordinates"]=np.column_stack((x,y))
                chr_db[variant]["var_info"]=np.array([v[2] for v in hits])
                chr_db[variant]["index"]=np.array([v[3] for v in hits])

                if args.DBSCAN:
                    db = DBSCAN(eps=args.epsilon, min_samples=args.min_pts).fit(chr_db[variant]["coordinates"])
                else:
                    db = DBSCAN(eps=args.bnd_distance, min_samples=2).fit(chr_db[variant]["coordinates"])

                core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                core_samples_mask[db.core_sample_indices_] = True
                labels = db.labels_
                #print variant
                n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

                unique_labels = set(labels)
                stats=[]
                #print the unique variants
                unique_xy=chr_db[variant]["coordinates"][  labels == -1 ]
                unique_index=chr_db[variant]["index"][  labels == -1 ]
    
                for j in range(0,len( chr_db[variant]["coordinates"][ labels == -1] ) ):
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
                                        
                #print the clusters
                for k in unique_labels:
                    class_member_mask = (labels == k)
                    xy = chr_db[variant]["coordinates"][class_member_mask & core_samples_mask]
                    indexes=chr_db[variant]["index"][class_member_mask & core_samples_mask]
                    if k == -1:
                        pass
                    elif args.DBSCAN:
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
                        variant_dictionary=fetch_index_variant(c,indexes)
                        chain_data=expand_chain(variant_dictionary,args.bnd_distance,args.overlap,args.ci)
                        clusters=cluster_variants(variant_dictionary,set(indexes),chain_data)
                        for clustered_variants in clusters:
                            f.write( vcf_line(clustered_variants,"cluster_{}".format(i),sample_IDs )+"\n")
                            i += 1
                        
    f.close()




 
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
    f.write( db_header()+"\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(sample_IDs)))
    f.close()
    export(args,sample_IDs)
