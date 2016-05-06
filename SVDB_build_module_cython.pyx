import SVDB_overlap_module
import random
import readVCF
import SVDB_merge_vcf_module_cython
#these modules are compiled using cython prior to constructing the db

#search through all variants and divide them into clusters
def generate_cluster(n,indexes,unclustered,chrA,chrB,args):
    changed=True
    cluster=set([n])
    tested=set([])
    
    #expand the cluster of variant n until the cluster cannot grow any more
    while changed:
        cluster_size=len(cluster)
        untested=cluster-tested
        for variant in untested:
            for i in indexes:
                hit=False
                #check if the variant is similar to any variant in the cluster, then add it
                if not i in cluster:
                    if not args.ci:
                        if chrA == chrB:
                            hit=SVDB_overlap_module.isSameVariation(unclustered[variant][0],unclustered[variant][1],unclustered[i][0],unclustered[i][1],args.overlap)
                        else:
                            hit=SVDB_overlap_module.precise_overlap(unclustered[variant][0],unclustered[variant][1],unclustered[i][0],unclustered[i][1],args.bnd_distance)
                    else:
                        hit=SVDB_overlap_module.ci_overlap_two_sided(unclustered[variant][0],unclustered[variant][1],unclustered[variant][2],unclustered[variant][3],unclustered[i][0],unclustered[i][1],unclustered[i][2],unclustered[i][3])
                    #print hit
                    if hit:
                        cluster=cluster | set([i])
        tested = tested | untested             
        if cluster_size == len(cluster):
            changed=False
    return(cluster)
    
#evaluate clusters and break large complex clusters into simple overlaped variants
def evaluate_cluster(cluster,samples,chrA,chrB,args):
    processed_clusters=[]
    cluster_list=list(cluster)
    total_list=cluster_list            
    while cluster_list:
        variant_index=random.randint(0,len(cluster_list)-1)
        variant=cluster_list[variant_index]
        del cluster_list[variant_index];
        overlapping=[]
        for j in total_list:
            hit = False
            if j != variant:
                if not args.ci:
                    if chrA == chrB:
                        hit=SVDB_overlap_module.isSameVariation(samples[variant][0],samples[variant][1],samples[j][0],samples[j][1],args.overlap)
                    else:
                        hit=SVDB_overlap_module.precise_overlap(samples[variant][0],samples[variant][1],samples[j][0],samples[j][1],args.bnd_distance)
                else:
                    hit=SVDB_overlap_module.ci_overlap_two_sided(samples[variant][0],samples[variant][1],samples[variant][2],samples[variant][3],samples[j][0],samples[j][1],samples[j][2],samples[j][3])
                if hit:
                    overlapping.append(j)
                    
        #delete all the matched variants, add a new entry into the processed cluster list
        overlapping=set(overlapping)
        cluster_list=list(set(cluster_list)-overlapping)
        processed_clusters.append(samples[variant]+[len(overlapping)+1])
        #add sample id of all the matched variants
        for j in list(overlapping):
            processed_clusters[-1][-2] += "|" + samples[j][-1]
      
    return(processed_clusters)

#get all the variants, and store them as a dictionary    
def get_variant_files(samples):
    chromosome_order=[]
    variant_dictionary={}
    Nsamples=0
    for sample in samples:
        Nsamples +=1
        sample_id=sample.split("/")[-1]
        sample_id=sample_id.replace(".vcf","")
        with open(sample) as variants:
            for line in variants:
                if(line[0] != "#"):
                    chrA,posA,chrB,posB,event_type =readVCF.readVCFLine(line);
                    if not chrA in variant_dictionary:
                        variant_dictionary[chrA]={}
                        if not chrA in chromosome_order:
                            chromosome_order.append(chrA)
                    if not chrB in variant_dictionary[chrA]:
                        variant_dictionary[chrA][chrB]={}
                    if not event_type in variant_dictionary[chrA][chrB]:
                        variant_dictionary[chrA][chrB][event_type]={}
                    if not sample_id in variant_dictionary[chrA][chrB][event_type]:
                        variant_dictionary[chrA][chrB][event_type][sample_id] =[]
                    CIA,CIB=SVDB_merge_vcf_module_cython.get_CIPOS_CEND(line.strip())
                    variant_dictionary[chrA][chrB][event_type][sample_id].append([int(posA),int(posB),CIA,CIB])
    return(variant_dictionary,chromosome_order,Nsamples)

