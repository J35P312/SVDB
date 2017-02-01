from __future__ import absolute_import
import sys, os, glob
from . import readVCF
from . import overlap_module
from operator import itemgetter
from . import merge_vcf_module_cython

def update_tag_data(tag_data,entry):
    if not entry["tag_id"] in tag_data:
        tag_data[entry["tag_id"]]=""
        tag_data[entry["tag_id"]] += "{}:{}-{}".format(entry["info"],entry["start"],entry["end"])
    else:
        tag_data[entry["tag_id"]] += "|{}:{}-{}".format(entry["info"],entry["start"],entry["end"])
    return tag_data


def main(args):

    bed_data={}
    tag_ids=[]
    
    for bed in args.file:
        tag_id=bed.split("/")[-1].replace(".txt","").replace(".bed","").replace(".tab","").replace(".csv","")
        tag_ids.append(tag_id)
        for line in open(bed):
            if line[0] == "#":
                continue
            content=line.split("\t")
            
            if not content[0].replace("chr","").replace("Chr","").replace("CHR","") in bed_data:
                bed_data[content[0].replace("chr","").replace("Chr","").replace("CHR","")] = []
        
            if len(content) > 3:
                entry={"chr":content[0].replace("chr","").replace("Chr","").replace("CHR",""),"start":int(content[1]),"end":int(content[2]),"info":content[3].strip(),"tag_id":tag_id}
            else:
                entry={"chr":content[0].replace("chr","").replace("Chr","").replace("CHR",""),"start":int(content[1]),"end":int(content[2]),"info":args.tag,"tag_id":tag_id}
                if not args.tag:
                    print ("error: the fourth column of the bed file was not found, use the --tag argument to set an id tag")
                    quit()
            bed_data[content[0].replace("chr","").replace("Chr","").replace("CHR","")].append(entry)

    

    for line in open(args.vcf):
        if line[0] == "#" and line[1] == "#":
            print line.strip()
            continue
        if line[0] == "#":
            for tag_id in tag_ids:
                print "##INFO=<ID={},Number=1,Type=String,Description=\"Bed annotation\">".format(tag_id)
            print line.strip()
            continue
        
        chrA,posA,chrB,posB,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
        if not chrA in bed_data or not chrB in bed_data:
            continue
            
        tag_data={}
        for entry in bed_data[chrA]:
            if entry["chr"] == chrA:
                if chrA != chrB:
                    if abs( posA- entry["start"]) <= args.bnd_distance or abs(posA-entry["end"]) <= args.bnd_distance:
                        tag_data=update_tag_data(tag_data,entry)
                        
                    elif( entry["start"] <= posA and posA <= entry["end"]):
                        tag_data=update_tag_data(tag_data,entry)
                else:
                    variant_length=posB-posA
                    entry_length=entry["end"] -entry["start"]
                    
                    start_end_distance=entry["end"]-posA
                    end_end_distance=entry["end"]-posB
                    
                    if start_end_distance < 0:
                        continue
                    elif end_end_distance > variant_length:
                        continue
                        
                    #the entire vcf variant fits within the bed entry
                    intersect= variant_length
                    #the entire bed entry fits withing the vcf entry
                    if(end_end_distance < 0 and start_end_distance > entry_length):
                        intersect= entry_length
                    #the vcf entry overlap the bed entry, but not fully, the end of the vcf > end of bed
                    elif(end_end_distance < 0):
                        intersect = variant_length +  end_end_distance
                    elif(start_end_distance > entry_length):
                        intersect = entry_length-end_end_distance
                    
                    
                    overlap=intersect/float(variant_length)
                    if overlap >= args.percentage:
                        tag_data=update_tag_data(tag_data,entry)
            elif entry["chr"] == chrB:
                if abs( posB- entry["start"]) <= args.bnd_distance or abs(posB-entry["end"]) <= args.bnd_distance:
                    tag_data=update_tag_data(tag_data,entry)
                elif( entry["start"] <= posB and posB <= entry["end"]): 
                    tag_data=update_tag_data(tag_data,entry)
        content=line.split("\t")
        if tag_data:
            for tag in tag_data:
                content[7]+= ";{}={}".format(tag,tag_data[tag])
        print("\t".join(content).strip() ) 
