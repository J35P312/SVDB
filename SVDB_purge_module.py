import SVDB_overlap_module
import readVCF

#delete specified samples from the db
def purge_sample(args):
    for line in open(args.db):
        if not "#" == line[0]:
            content=line.strip().split("\t")
            info=content[7]
            samples=info.split(";SAMPLES=")[-1];
            samples=samples.split("|");

            sample_set=set(samples)-set(args.samples)
            hits=len(sample_set)
            db_size=info.split("NSAMPLES=")[-1]
            db_size=float(db_size.split(";")[0])-len(args.samples)
            frequency=hits/db_size
            samples="|".join(list(sample_set))
            variant_info="OCC={};NSAMPLES={};FRQ={};SAMPLES={}".format(hits,db_size,frequency,samples)
            info=info.split("OCC")[0];
            info+= variant_info;

            content[7]=info
            if hits:
                print("\t".join(content))

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
