import SVDB_overlap_module

def retrieve_key(line,key):
    item = False
    if key + "=" in line:
        item=line.strip().split(key+"=")[-1].split(";")[0]
        if(len(item) == len(line.strip())):
            return(False)
    return(item)

def get_CIPOS_CEND(query_variant):
    ciA_query=[0,0]
    CIPOS=retrieve_key(query_variant[-1],"CIPOS")
    if CIPOS:
        CIPOS=CIPOS.split(",")
        if len(CIPOS) == 2:
            ciA_query=[int(CIPOS[0]),int(CIPOS[1])]
        else:
            ciA_query=[int(CIPOS[0]),int(CIPOS[0])]

    ciB_query=[0,0]
    CIPOS=retrieve_key(query_variant[-1],"CIEND")
    if CIPOS:
        CIPOS=CIPOS.split(",")
        if len(CIPOS) == 2:
            ciB_query=[int(CIPOS[0]),int(CIPOS[1])]
        else:
            ciB_query=[int(CIPOS[0]),int(CIPOS[0])]

    return(ciA_query,ciB_query)

def find_ci(query_variant,db_variant):

    ciA_query,ciB_query=get_CIPOS_CEND(query_variant)
    ciA_db,ciB_db=get_CIPOS_CEND(db_variant)

    return(ciA_query,ciB_query,ciA_db,ciB_db)

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

def merge(variants,ci,overlap,bnd_distance,no_intra):
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
                    if not no_intra or variants[chrA][i][-3] != variants[chrA][j][-3]:
                        if not ci:
                            overlap=SVDB_overlap_module.variant_overlap(chrA,variants[chrA][i][0],variants[chrA][i][2],variants[chrA][i][3],variants[chrA][j][2],variants[chrA][j][3],overlap,bnd_distance)
                        else:
                            ciA_query,ciB_query,ciA_db,ciB_db=find_ci(variants[chrA][i],variants[chrA][j])
                            overlap=SVDB_overlap_module.ci_overlap(variants[chrA][i][2],variants[chrA][i][3],ciA_query,ciB_query,variants[chrA][j][2],variants[chrA][j][3],ciA_db,ciB_db)
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
            if not line[0] in to_be_printed:
                to_be_printed[line[0]]=[]
            to_be_printed[line[0]].append(line)
            i +=1

    return(to_be_printed)
