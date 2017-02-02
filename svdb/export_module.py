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
    A='SELECT var ,chrA , chrB ,posA ,ci_A_lower ,ci_A_upper ,posB ,ci_B_lower ,ci_B_upper FROM SVDB WHERE idx == \'{}\' '.format(index)
    hits = c.execute(A)
    variant={}
    for hit in hits:
        variant["type"]=hit[0]
        variant["chrA"]=hit[1]
        variant["chrB"]=hit[2]
        variant["posA"]=int( hit[3] )
        variant["ci_A_start"]= int(hit[4])
        variant["ci_A_end"]= int(hit[5])
        variant["posB"]=int( hit[6] )
        variant["ci_B_start"]= int(hit[7])
        variant["ci_B_end"]= int(hit[8])
        
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
    headerString+="##INFO=<ID=NSAMPLES,Number=1,Type=Int,Description=\"the number of samples within the database\">\n";
    headerString+="##INFO=<ID=VARIANTS,Number=1,Type=Int,Description=\"a| separated list of the positions of the clustered variants\">\n";
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

def expand_chain(chain,c,distance,overlap,ci):
    tested=set([])
    chain_data={}
    while not tested == chain:
        to_be_tested=chain-tested
        for index in to_be_tested:
            chain_data[index]=[]
            variant=fetch_index_variant(c,index)
            selection ="idx"
            if variant["chrA"] == variant["chrB"]:
                selection = "posA, posB, idx"
            if not ci:
                A='SELECT {} FROM SVDB WHERE var == \'{}\' AND chrA == \'{}\' AND chrB == \'{}\' AND posA <= {} AND posA >= {} AND posB <= {} AND posB >= {}'.format(selection,variant["type"],variant["chrA"],variant["chrB"],variant["posA"]+distance, variant["posA"] -distance,variant["posB"] + distance, variant["posB"]-distance)
                #print A
            else:
                A='SELECT idx FROM SVDB WHERE var == \'{}\' AND chrA == \'{}\' AND chrB == \'{}\' AND posA <= {} AND posA >= {} AND posB <= {} AND posB >= {}'.format(variant["type"],variant["chrA"],variant["chrB"],variant["posA"]+variant["ci_A_end"], variant["posA"] -variant["ci_A_start"],variant["posB"] + variant["ci_B_end"], variant["posB"] - variant["ci_B_start"])
            
            hits = c.execute(A)
            for hit in hits:
                if not ci and variant["chrA"] == variant["chrB"]:
                    var={}
                    var["posA"]=int( hit[0] )
                    var["posB"]=int( hit[1] )
                    var["index"]=int( hit[2] )
                    similar=overlap_module.isSameVariation(variant["posA"],variant["posB"],var["posA"],var["posB"],overlap,distance)
                    if similar:
                        chain_data[index].append(var["index"])
                        chain.add( int( var["index"]) )
                    
                else:
                    chain_data[index].append(int( hit[0] ))
                    chain.add( int( hit[0] ) )

            tested.add(index)             
            #print"{} {} {} {} {}".format(variant["type"],variant["chrA"],variant["chrB"],variant["posA"],variant["posB"])
    return([chain,chain_data])

def cluster(c,chain,chain_data,distance,overlap,ci):
    clusters=[]
    for index in sorted(chain_data, key=lambda index: len(chain_data[index]), reverse=True):
        if not chain:
            break
        if not index in chain:
            continue
        for var in chain_data[index]:
            if var in chain:
                chain.remove(var)
        variant_dictionary=fetch_cluster_variant(c,chain_data[index])
        variant=fetch_index_variant(c,index)
       
        
        
        
        clusters.append( [variant,variant_dictionary] )      
    return(clusters)

def generate_chains(db,prefix,distance,overlap,ci,sample_IDs,args):

    
    f = open(prefix+".vcf",'a')
    conn = sqlite3.connect(db)
    if args.memory:
        memory_db=sqlite3.connect(':memory:')
        db_dump="".join(line for line in conn.iterdump())
        memory_db.executescript(db_dump)
        conn.close
        c = memory_db.cursor()
    else:
        c=conn.cursor()

    chains=[]
    A="SELECT MAX(idx) FROM SVDB"
    number_of_variants=0
    for hit in c.execute(A):
        number_of_variants=int(hit[0])

    variant_list=set( range(0,number_of_variants) )
    summed=0;
    i =1
    while variant_list:
        chain = set([ variant_list.pop() ])
        #find similar variants
        chain, chain_data= expand_chain(chain,c,distance,overlap,ci)
        for index in chain:
            if index in variant_list:
                variant_list.remove(index)

        clusters=cluster(c,chain,chain_data,distance,overlap,ci)
        for clustered_variants in clusters:        
            f.write( vcf_line(clustered_variants,"cluster_{}".format(i),sample_IDs )+"\n")
            i += 1
            
            
    f.close()
def dbscan_export(args,sample_IDs):

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
    i=0;
    for chrA in chrA_list:
        for chrB in chrB_list:
            chr_db={}            
            for hit in c.execute('SELECT posA,posB,sample,idx,var FROM SVDB WHERE chrA == \'{}\' AND chrB == \'{}\''.format(chrA,chrB)):
                if not hit[-1] in chr_db:
                    chr_db[ hit[-1] ]={}
                    chr_db[ hit[-1] ]["coordinates"]=[]
                    chr_db[ hit[-1] ]["var_info"]=[]
                    chr_db[ hit[-1] ]["index"]=[]
                
                chr_db[ hit[-1] ]["coordinates"].append([hit[0],hit[1]])
                chr_db[ hit[-1] ]["var_info"].append(hit[2])
                chr_db[hit[-1]]["index"].append(hit[-2])
                  
            for variant in chr_db:
                chr_db[variant]["coordinates"]=np.array(chr_db[variant]["coordinates"])
                chr_db[variant]["var_info"]=np.array(chr_db[variant]["var_info"])
                chr_db[variant]["index"]=np.array(chr_db[variant]["index"])
                  
                db = DBSCAN(eps=500, min_samples=2).fit(chr_db[variant]["coordinates"])
                
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
                    else:
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
    if args.DBSCAN:
        dbscan_export(args,sample_IDs)
    else:
        generate_chains(args.db,args.prefix,args.bnd_distance,args.overlap,args.ci,sample_IDs,args )
