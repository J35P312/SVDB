import SVDB_overlap_module
import SVDB_build_module
import SVDB_query_module
import glob,os
import random
import numpy

def build_db(args,samples,prefix):
    args.files=samples
    args.folder=False
    args.overlap=0.9
    args.bnd_distance=500
    args.prefix=prefix;
    SVDB_build_module.main(args)

def query_db(args,sample,db,prefix):
    args.prefix=prefix;
    args.invert=False
    args.query_vcf=sample
    args.db=db
    args.hit_tag="OCC"
    args.frequency_tag="FRQ"
    args.prefix=prefix
    args.no_var=False
    SVDB_query_module.main(args)

def get_frequencies(queried_vcf):
    frequency_histogram={}
    lines=0;
    with open(queried_vcf) as f:
        for line in f:

            content=line.split("\t")
            if not line[0] == "#":
                lines+=1
                info=line.split("\t")[7]
                frequency=info.split(";OCC=")[-1];
                frequency=float(frequency.split(";")[0])
                if not frequency in frequency_histogram:
                    frequency_histogram[frequency]=0
                frequency_histogram[frequency] += 1
    return(frequency_histogram,lines)

def clear_db(prefix):
    os.remove(prefix+".db.vcf")

def similarity_matrix(args,samples):
    first= True;
    #print the header
    for sample in sorted(samples):
        if first:
            first=False
            header=sample
        else:
            header += "," + sample
    print(header)

    #generate a db per sample
    for sample in samples:
        build_db(args,[samples[sample]],sample)

    #each sample will query each db
    for query_sample in samples:
        for db_sample in samples:
            if not os.path.exists(query_sample+"_"+db_sample+"_query.vcf"):
                query_db(args,samples[query_sample],db_sample+".db.vcf",query_sample+"_"+db_sample)

    #generate each row of the matrix
    a=range(len(samples))
    matrix= [[i*0 for i in a] for i in a]
    i=0;
    for row_sample in sorted(samples):
        row=[]
        for column_sample in sorted(samples):
            frequency_histogram,variants=get_frequencies(row_sample+"_"+column_sample+"_query.vcf")
            hits=0;
            for frequency in frequency_histogram:
                if not frequency == 0:
                     hits +=frequency_histogram[frequency]
            row.append(hits)
            #os.remove(row_sample+"_"+column_sample+"_query.vcf")
        norm=row[i]
        for j in range(0,len(row)):
            matrix[i][j]=round(row[j]/float(norm),2)
        i+=1;

    #calculate the mean value of the queries
    for i in range(0,len(matrix)):
        for j in range(0,len(matrix)):
            avg=(matrix[i][j]+matrix[j][i])/2
            matrix[i][j]=avg
            matrix[j][i]=avg

    i=0;
    for row_sample in sorted(samples):  
        row=[row_sample]
        for j in range(0,len(samples)):
            row.append(str(matrix[i][j]))
        print(",".join(row))
        i +=1

    #clear all query files
    for row_sample in sorted(samples):
        for column_sample in sorted(samples):
            if os.path.exists(row_sample+"_"+column_sample+"_query.vcf"):
                os.remove(row_sample+"_"+column_sample+"_query.vcf")
    #clear all the db files
    for sample in samples:
        clear_db(sample)

def sample_hist(args,samples):
    sample_list=[]
    for sample in samples:
        sample_list.append(samples[sample])
    print("k,avg_unique,std_unique")
    if not args.k:
        l=range(10,len(samples))
        args.k=l[0::10]
        if not len(samples) in args.k:
            args.k.append(len(samples))
    frequency_hist={}
    prefix="tmp"+"_generate_hist"

    for k in args.k:
        ones=[]
        hist={}
        for n in range(0,args.n):
            db_samples=random.sample(sample_list,k)
            
            build_db(args,db_samples,prefix)
            #use the database to query each sample
            for sample in db_samples:
                query_db(args,sample,prefix+".db.vcf","tmp")
                frequency_hist[sample],variants=get_frequencies("tmp_query.vcf")                
                os.remove("tmp_query.vcf")

                if 1 in frequency_hist[sample]:
                    ones.append(frequency_hist[sample][1]/float(variants))

                for val in frequency_hist[sample]:
                    if not val in hist:
                        hist[val] = []
                    hist[val].append(frequency_hist[sample][val]/float(variants))
                
        
        print("{},{},{}".format(k,numpy.average(ones),numpy.std(ones)))
        f=open("SVDB_hist_{}.csv".format(k),"w")
        f.write("frequency,variant_frequency\n")
        for val in sorted(hist):
            f.write("{},{}\n".format(int(val),numpy.mean(hist[val])))
    os.remove(prefix+".db.vcf")  
            

def main(args):
    if args.folder:
        samples=glob.glob("{}/*.vcf".format( os.path.abspath(args.folder) ) )
    elif args.files:
        samples=args.files
    sample_dictionary={}
    for sample in samples:
        name=sample.split("/")[-1]
        name=name.replace(".vcf" ,"")
        sample_dictionary[name]=sample

    if args.similarity_matrix:
        similarity_matrix(args,sample_dictionary);
    elif(args.sample_hist):
        sample_hist(args,sample_dictionary);
    else:
        print("choose the sample_hist or similarity matrix mode")
