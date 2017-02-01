from __future__ import absolute_import
from . import overlap_module
from . import readVCF

#delete specified samples from the db

def purge_samples_info_field(sample_info,samples):
    new_sample_fields=""
    for info in sample_info:
        if not info.split(":")[0] in samples:
            new_sample_fields += info + "|"

    return(new_sample_fields)
    


def purge_sample(args):
    db_size=0;
    sample_pos_vector=[]
    for line in open(args.db):
        if not "#" == line[0]:
            content=line.strip().split("\t")
            info=content[7]
            samples=info.split(";VARIANTS=")[-1];
            samples=samples.split("|");
            samples=purge_samples_info_field(samples,args.samples)

            sample_set=set(samples)-set(args.samples)
            
            
            
            for i in sorted(sample_pos_vector, reverse=True):
                del content[i]
               
            hits=content.count("0/0")
            hits=db_size-hits
            frequency=round( hits/float(db_size) , 2)
            variant_info="NSAMPLES={};OCC={};FRQ={};VARIANTS={}".format(int(db_size),hits,frequency,samples)
            info=info.split("NSAMPLES")[0];
            info+= variant_info;

            content[7]=info
            if hits:
                print("\t".join(content))

        else:
            if("#CHROM\t" in line):
                content=line.split("\t")
                #remove the position and info fields from the count
                for sample in args.samples:
                    if sample in content:
                        sample_pos_vector.append(content.index(sample))
                
                for i in sorted(sample_pos_vector, reverse=True):
                    del content[i]
                
                db_size=len(content)-9
                print( "\t".join(content).strip() )         
            else:
                print(line.strip())
#remove variants if they are found in the vcf
def purge_variants(args):
    variants={}
    for line in open(args.vcf):
        if not line[0] == "#":
            chrA,chrApos_query,chrB,chrBpos_query,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
            if not chrA in variants:
                variants[chrA]={}
            if not chrB in variants[chrA]:
                variants[chrA][chrB]={}
            if not event_type in variants[chrA][chrB]:
                variants[chrA][chrB][event_type] = []
            variants[chrA][chrB][event_type].append([chrApos_query,chrBpos_query])

    for db_line in open(args.db):
        if db_line[0] == "#":
            print(db_line.strip())
        else:

            line_found= False
            chrA,chrApos_db,chrB,chrBpos_db,event_type,INFO,FORMAT =readVCF.readVCFLine(db_line);
            if chrA in variants:
                if chrB in variants[chrA]:
                    if event_type in variants[chrA][chrB]:
                        for var in variants[chrA][chrB][event_type]:
                            overlap=SVDB_overlap_module.variant_overlap(chrA,chrB,var[0],var[1],chrApos_db,chrBpos_db,args.overlap,args.bnd_distance)
                            if overlap:
                                line_found=True

            if not line_found:
                print(db_line.strip())

def main(args):
    if args.samples or args.file:

        if args.file:
            args.samples=[]
            for line in open(args.file):
                if not line.strip() == "":
                    args.samples.append(line.strip())

        purge_sample(args)
    elif args.vcf:
        purge_variants(args)
