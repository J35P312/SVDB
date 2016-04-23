import readVCF
import SVDB_overlap_module

def merge_csq(info,csq):
    var_csq=info.split("CSQ=")[-1].split(";")[0]
    csq.append(var_csq);
    effects=set([])
    for csq_field in csq:
        effects= effects | set(csq_field.split(","))
    CSQ="CSQ="+",".join(list(effects))
    pre_CSQ=info.split("CSQ=")[0]
    post_CSQ=info.split("CSQ=")[-1].split(";")
    if len(post_CSQ) == 1:
        post_CSQ=""
    else:
        post_CSQ=";"+";".join(post_CSQ[1:])

    info=pre_CSQ+CSQ+post_CSQ

    return(info)

def main(args):
    variants={}
    #add all the variants into a dictionary
    i=0;
    for line in open(args.vcf):
        if line[0] == "#":
            if line[1] == "#":
                print(line.strip())
            else:
                print("##INFO=<ID=VARID,Number=1,Type=String,Description=\"The variant ID of merged samples\">")
                print(line.strip())
        else:
            chrA,posA,chrB,posB,event_type =readVCF.readVCFLine(line);
            if not chrA in variants:
                variants[chrA]=[]
            variants[chrA].append([chrB,event_type,posA,posB,i,line.strip()])
            i+=1;

    #search for similar variants
    to_be_printed={}
    for chrA in variants:
        i=0;

        while i < len(variants[chrA]):
            merge=[]
            csq=[]
            j=i+1;
            while j < len(variants[chrA]):
                    
                if variants[chrA][i][0] == variants[chrA][j][0] and variants[chrA][i][1] == variants[chrA][j][1]:
                    overlap=SVDB_overlap_module.variant_overlap(chrA,variants[chrA][i][0],variants[chrA][i][2],variants[chrA][i][3],variants[chrA][j][2],variants[chrA][j][3],args.overlap,args.bnd_distance)
                    if overlap:
                        #add similar variants to the merge list and remove them
                        merge.append(variants[chrA][j][-1].split("\t")[2])
                        if variants[chrA][i][0] != chrA and "CSQ=" in variants[chrA][j][-1]:
                            info=variants[chrA][j][-1].split("\t")[7]
                            csq.append(info.split("CSQ=")[-1].split(";")[0])
                        del variants[chrA][j]
                        j += -1
                j+=1
            line=variants[chrA][i][-1].split("\t")
            line[7] += ";VARID=" + "|".join(merge)
            if csq:
                line[7]=merge_csq(line[7],csq)
            to_be_printed[variants[chrA][i][-2]]="\t".join(line)
            i +=1
    #print the variants in similar order as the input
    for var in sorted(to_be_printed):
        print(to_be_printed[var])
