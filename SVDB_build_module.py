import sys, os, glob
import argparse
import readVCF
import SVDB_overlap_module
import subprocess
import sqlite3

def biggest_cluster(chain,chain_data):
    biggest=-1
    biggest_cluster=-1
    for index in chain:
        if len(chain_data[index]) > biggest:
            biggest=len(chain_data[index])
            biggest_cluster=index
    return(biggest_cluster)

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


def print_to_vcf(cluster,id_tag,sample_IDs):
    info_field="SVTYPE={};".format(cluster[0]["type"])
    vcf_line=[]
    vcf_line.append( cluster[0]["chrA"] )
    vcf_line.append( str(cluster[0]["posA"]) )
    vcf_line.append( id_tag )
    vcf_line.append( "N" )
    if cluster[0]["chrA"] == cluster[0]["chrB"]:
        vcf_line.append( cluster[0]["type"] )
        info_field += "END={};SVLEN={};".format(cluster[0]["posB"],abs(cluster[0]["posA"]-cluster[0]["posB"]))
    else:
        vcf_line.append( "N[{}:{}[".format(cluster[0]["chrB"],cluster[0]["posB"]) )
        
        
    sample_set=set([])
    for variant in cluster[1]:
        sample_set = sample_set | set([cluster[1][variant]["sample_id"]])
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
    for sample in zygosity_list:
        format.append(zygosity_list[sample])
    vcf_line.append("GT")
    vcf_line.append("\t".join(format))
    return( "\t".join(vcf_line) )

def populate_db(args):
    sample_IDs=[]
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    tableListQuery = "SELECT name FROM sqlite_master WHERE type=\'table\'"
    c.execute(tableListQuery)
    tables = map(lambda t: t[0], c.fetchall())
    new_db = False
    if not "SVDB" in tables:
        new_db = True
        A="CREATE TABLE SVDB (var TEXT,chrA TEXT, chrB TEXT,posA INT,posB INT, sample TEXT, idx INT)"
        c.execute(A)
    else:
        A="DROP INDEX SV "
        c.execute(A)
        A="DROP INDEX IDX"
        c.execute(A)
         
    idx=0;
    #populate the tables
    for vcf in args.files:
        sample_name=vcf.split("/")[-1]
        sample_name=sample_name.replace(".vcf","")
        sample_name=sample_name.replace(".","_")
        sample_IDs.append(sample_name)
        A='SELECT sample FROM SVDB WHERE sample == \'{}\' '.format(sample_name)
        hits=[]
        for hit in c.execute(A):
            hits.append(hit)
            
        if hits:
            continue        
        
        var =[]
        for line in open(vcf):
            if line.startswith("#"):
                continue
            chrA,posA,chrB,posB,event_type,INFO,FORMAT = readVCF.readVCFLine(line)
            var.append((event_type,chrA,chrB,posA,posB,sample_name,idx))
            idx += 1;
            #insert EVERYTHING into the database, the user may then query it in different ways(at least until the DB gets to large to function properly)
        c.executemany('INSERT INTO SVDB VALUES (?,?,?,?,?,?,?)',var)     
        

    A="CREATE INDEX SV ON SVDB (var,chrA, chrB, posA , posA,posB,posB)"
    c.execute(A)
    A="CREATE INDEX IDX ON SVDB (idx)"
    c.execute(A)
    #test query
    #A = "EXPLAIN QUERY PLAN SELECT * FROM SVDB WHERE var == \'DEL\' AND chrA == \'1\' AND chrB == \'1\' AND posA <= 1 AND posA >= 1 AND posB <= 1 AND posB >= 1"
    #print "query plan:"
    #for hit in c.execute(A):
    #    print hit
        
    conn.commit()
    conn.close()
    return(sample_IDs)

def expand_chain(chain,c,distance,overlap,ci):
    tested=set([])
    chain_data={}
    while not tested == chain:
        
        to_be_tested=chain-tested
        for index in to_be_tested:
            chain_data[index]={}
            A='SELECT * FROM SVDB WHERE idx == \'{}\' '.format(index)
            hits = c.execute(A)
            variant={}
            for hit in hits:
                variant["type"]=hit[0]
                variant["chrA"]=hit[1]
                variant["chrB"]=hit[2]
                variant["posA"]=int( hit[3] )
                variant["ci_A_start"]=0
                variant["ci_A_end"]=0
                variant["posB"]=int( hit[4] )
                variant["ci_B_start"]=0
                variant["ci_B_end"]=0
                
            if not ci:
                A='SELECT * FROM SVDB WHERE var == \'{}\' AND chrA == \'{}\' AND chrB == \'{}\' AND posA <= {} AND posA >= {} AND posB <= {} AND posB >= {}'.format(variant["type"],variant["chrA"],variant["chrB"],variant["posA"]+distance, variant["posA"] -distance,variant["posB"] + distance, variant["posB"]-distance)
            else:
                A='SELECT * FROM SVDB WHERE var == \'{}\' AND chrA == \'{}\' AND chrB == \'{}\' AND posA <= {} AND posA >= {} AND posB <= {} AND posB >= {}'.format(variant["type"],variant["chrA"],variant["chrB"],variant["posA"]+variant["ci_A_end"], variant["posA"] -variant["ci_A_start"],variant["posB"] + variant["ci_B_start"], variant["posA"]-variant["ci_B_end"])
            
            similar_variants=[]
            hits = c.execute(A)
            for hit in hits:
                var={}
                var["type"]=hit[0]
                var["chrA"]=hit[1]
                var["chrB"]=hit[2]
                var["posA"]=int( hit[3] )
                var["ci_A_start"]=0
                var["ci_A_end"]=0
                var["posB"]=int( hit[4] )
                var["ci_B_start"]=0
                var["ci_B_end"]=0
                var["sample_id"]=hit[-2]
                var["index"]=int( hit[-1] )
                similar_variants.append(var)
            
            #ANCIENT CODE!
            
            for sv in similar_variants:
                if not args.ci:
                    if sv["chrA"] == sv["chrB"]:
                        hit=SVDB_overlap_module.isSameVariation(variant["posA"],variant["posB"],sv["posA"],sv["posB"],overlap,distance)
                    else:
                        hit=SVDB_overlap_module.precise_overlap(variant["posA"],variant["posB"],sv["posA"],sv["posB"],distance)
                else:
                    hit=SVDB_overlap_module.ci_overlap_two_sided(variant["posA"],variant["posB"],[variant["ci_A_start"],variant["ci_A_end"]],[variant["ci_B_start"],variant["ci_B_end"]],sv["posA"],sv["posB"],[ sv["ci_A_start"],sv["ci_A_end"] ],[sv["ci_B_start"],sv["ci_B_end"]])
                    
                        #print hit
                if hit:
                    chain_data[index][sv["index"]]=sv
                    chain=chain | set([ sv["index"] ])
            tested = tested | set([index])             
            #print"{} {} {} {} {}".format(variant["type"],variant["chrA"],variant["chrB"],variant["posA"],variant["posB"])
    return([chain,chain_data])

def cluster(c,chain,chain_data,distance,overlap,ci):
    clusters=[]
    A=chain  
    while A:
        var_A=biggest_cluster(A,chain_data)
        current_cluster=[var_A]
        variant=chain_data[var_A][var_A]
        selected=[var_A]
        for var in chain_data[var_A]:
            selected.append(var)
     
        clusters.append( [ variant,chain_data[var_A]])
        A=A-set(selected)
        
    
    return(clusters)

def generate_chains(db,distance,overlap,ci,sample_IDs):

    #memory_db=sqlite3.connect(':memory:')
    conn = sqlite3.connect(db)
    #db_dump="".join(line for line in conn.iterdump())
    #memory_db.executescript(db_dump)
    #conn.close
    c = conn.cursor()
    #c = memory_db.cursor()

    chains=[]
    A="SELECT MAX(idx) FROM SVDB"
    number_of_variants=0
    for hit in c.execute(A):
        number_of_variants=int(hit[0])

    variant_list=set( range(0,number_of_variants) )
    summed=0;
    i =1
    while variant_list:
        chain= set([ min(variant_list) ])
        #find similar variants
        chain, chain_data= expand_chain(chain,c,distance,overlap,ci)
        clusters=cluster(c,chain,chain_data,distance,overlap,ci)
        for clustered_variants in clusters:
            print print_to_vcf(clustered_variants,"cluster_{}".format(i),sample_IDs )
            i += 1
        variant_list = variant_list - chain
        #chains.append(chain)
        #summed+=len(chains[-1])
        #print "{} {} {} {}".format(len(chains), len(chains[-1]),summed,number_of_variants) 
        
    return(chains)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("""This scripts takes as input any number of vcf files and generates a structural variant database""")
    parser.add_argument('--db', help="output database", type=str,default="SVDB.db")        
    parser.add_argument('--build'       , help="create a db", required=False, action="store_true")
    parser.add_argument('--no_merge'       , help="skip the merging of variants", required=False, action="store_true")
    parser.add_argument('--ci', help="overides overlap and bnd_distance,determine hits based on the confidence interval of the position fo the variants(0 if no CIPOS or CIEND is vailable)", required=False, action="store_true")
    parser.add_argument('--bnd_distance', type=int,default= 2500,help="the maximum distance between two similar precise breakpoints(default = 2500)")
    parser.add_argument('--overlap', type=float, default = 0.8,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.8")
    parser.add_argument('--files' , type=str, nargs='*', help="create a db using the specified vcf files(cannot be used with --folder)")
    parser.add_argument('--folder', type=str, help="create a db using all the vcf files in the folders")
    parser.add_argument('--prefix', type=str,default=None ,help="the prefix of the output file, default = print to stdout")
    
    args = parser.parse_args()
    if(args.files):
        sample_IDs = populate_db(args)
    elif(args.folder):
        vcf_folder = glob.glob(os.path.join(args.folder,"*.vcf"));
        test=[];
        args.files=vcf_folder
        sample_IDs = populate_db(args)
    else:
        print("error: use the vcf or folder option to select the input vcf")
        sys.exit()
    print( db_header() )
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format("\t".join(sample_IDs)))
    chains=generate_chains(args.db,args.bnd_distance,args.overlap,args.ci,sample_IDs )
