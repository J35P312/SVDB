import sys, os, glob
import SVDB_overlap_module
import SVDB_build_module_cython


def clustered(variant_dictionary,chrA,chrB,args):
    clustered_variants=[]
    unclustered=[]
    for sample in variant_dictionary:
        for variant in variant_dictionary[sample]:
            unclustered.append(variant+[sample])
    #divide the variants into clusters until all the variants belong to a cluster
    clusters=[]
    indexes=set(range(0,len(unclustered)))
    while len(indexes):
        new_cluster=SVDB_build_module_cython.generate_cluster(min(indexes),indexes,unclustered,chrA,chrB,args)
        clusters.append(new_cluster)
        #delete all the indexes belonging to the last cluster
        indexes=indexes - new_cluster
    #now merge the variants in the clusters
    for cluster in clusters:
        #no merging if there is only one variant in the cluster, if there are only 2 variants, pick represent the overlaping variants as the one having the lowest index
        if len(cluster) < 3:
            clustered_variants.append(unclustered[min(cluster)]+[len(cluster)])
            if len(cluster) > 1:
            	clustered_variants[-1][-2] +="|" + unclustered[max(cluster)][-1]
        else:
            clustered_variants += SVDB_build_module_cython.evaluate_cluster(cluster,unclustered,chrA,chrB,args)
    
    return(clustered_variants);

#the vcf header of the db
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
    headerString+="##INFO=<ID=SAMPLES,Number=1,Type=Int,Description=\"the ID of the samples carrying this variant\">\n";
    headerString+="##INFO=<ID=FRQ,Number=1,Type=Float,Description=\"the frequency of the vriant\">\n";
    headerString+="##INFO=<ID=samples,Number=1,Type=String,Description=\"the sample ID of each sample having this variant, seprated using |\">\n";
    headerString+="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    return(headerString)

def main(args):
    print(db_header())
    #store all the variants in a dictionary
    if args.build:
        samples=glob.glob("{}/*.vcf".format( os.path.abspath(args.build) ) )
    elif args.files:
        samples=args.files
                
    variant_dictionary,chromosome_order,Nsamples=SVDB_build_module_cython.get_variant_files(samples)
    variant_dictionary=SVDB_build_module_cython.clear_duplicate(variant_dictionary,args)                
    processed_variants={}
    #compute the frequency of the variants
    for chromosomeA in chromosome_order:
        processed_variants[chromosomeA]={}
        for chromosomeB in variant_dictionary[chromosomeA]:
            processed_variants[chromosomeA][chromosomeB]={}
            for event_type in variant_dictionary[chromosomeA][chromosomeB]:
                processed_variants[chromosomeA][chromosomeB][event_type]=clustered(variant_dictionary[chromosomeA][chromosomeB][event_type],chromosomeA,chromosomeB,args)
    #print the variants in the same order as the contig order in the first vcf file
    ID=0 
    for chromosomeA in chromosome_order:
        for chromosomeB in processed_variants[chromosomeA]:
            for variant_type in processed_variants[chromosomeA][chromosomeB]:
                for variant in processed_variants[chromosomeA][chromosomeB][variant_type]:
                    ID+=1
                    variant_tag=""
                    #convert TDUP etc into DUP
                    posA=variant[0]
                    posB=variant[1]
                    INFO="SVTYPE={};OCC={};NSAMPLES={};FRQ={};SAMPLES={}".format(variant_type,variant[-1],Nsamples,variant[-1]/float(Nsamples),variant[-2])
                    if(variant_type == "BND"):
                        variant_tag="N[{}:{}[".format(chromosomeB,posB)
                    else:
                        variant_tag= "<" + variant_type + ">"                        
                        INFO= "END={};".format(posB)+INFO
                        
                    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chromosomeA,posA,ID,"N",variant_tag,".","PASS",INFO))    
