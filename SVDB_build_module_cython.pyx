import SVDB_overlap_module
import random
import readVCF
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
                    #print unclustered[variant]
                    #print unclustered[i]
                    if chrA == chrB:
                        #print "test"
                        #print args.overlap
                        hit=SVDB_overlap_module.isSameVariation(unclustered[variant][0],unclustered[variant][1],unclustered[i][0],unclustered[i][1],args.overlap)
                    else:
                        hit=SVDB_overlap_module.precise_overlap(unclustered[variant][0],unclustered[variant][1],unclustered[i][0],unclustered[i][1],args.bnd_distance)
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
                if chrA == chrB:
                    hit=SVDB_overlap_module.isSameVariation(samples[variant][0],samples[variant][1],samples[j][0],samples[j][1],args.overlap)
                else:
                    hit=SVDB_overlap_module.precise_overlap(samples[variant][0],samples[variant][1],samples[j][0],samples[j][1],args.bnd_distance)
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
#clear duplicate entries from samples, for example when multiple callers have found the same variant 
def clear_duplicate(variant_dictionary,args):
       
    for chromosomeA in variant_dictionary:
        for chromosomeB in variant_dictionary[chromosomeA]:
            for event_type in variant_dictionary[chromosomeA][chromosomeB]:
                for sample_id in variant_dictionary[chromosomeA][chromosomeB][event_type]:
                    sample_variants=variant_dictionary[chromosomeA][chromosomeB][event_type][sample_id]
                    i=0;
                    while i < len(sample_variants):
                        j=0;
                        while j < len(sample_variants):
                            if j != i:
                                hit = False
                                if chromosomeA == chromosomeB:
                                        #if sample_variants[i][1] >= sample_variants[j][0] and sample_variants[j][1] >= sample_variants[i][0]:
                                        hit=SVDB_overlap_module.isSameVariation(sample_variants[i][0],sample_variants[i][1],sample_variants[j][0],sample_variants[j][1],0.9)
                                else:
                                    hit=SVDB_overlap_module.precise_overlap(sample_variants[i][0],sample_variants[i][1],sample_variants[j][0],sample_variants[j][1],500)
                                if hit:
                                    del sample_variants[j]
                                    j += -1
                            j+=1 
                        i+=1
                    variant_dictionary[chromosomeA][chromosomeB][event_type][sample_id]=sample_variants
    return(variant_dictionary)

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
                    variant_dictionary[chrA][chrB][event_type][sample_id].append([int(posA),int(posB)])
    return(variant_dictionary,chromosome_order,Nsamples)

