from __future__ import absolute_import
from . import readVCF
from . import merge_vcf_module_cython
import sys

def print_header(vcf_list,vcf_dictionary,args,command_line):
    header={};
    header["ALT"]={}
    header["INFO"]={}
    header["FILTER"]={}
    header["FORMAT"]={}
    header["CONTIGS"]=[]
    contigs=False
    reference=""
    subheader={}
    contigs_list=[]
    columns=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    sample_order={}
    print("##fileformat=VCFv4.1")
    print("##source=MergeVCF")
    print "##SVDB_version={} cmd=\"{}\"".format(args.version," ".join(sys.argv))
    samples=[]
    first=True
    first_vcf_header=""
    for vcf in vcf_list:
        for line in open(vcf):
            if(line[0] == "#"):
                if("#CHROM\tPOS" in line):
                    if first:
                        first_vcf_header=line.strip()
                        first=False
                    vcf_columns=line.strip().split("\t")
                    for column in vcf_columns:
                        if not column in columns:
                            
                            columns.append(column)                    
                            if not column in sample_order and column != "FORMAT":
                                sample_order[column]={}
                                samples.append(column)

                    if len(vcf_columns) > 8 and not args.same_order:
                        i=0
                        for sample in vcf_columns[9:]:                           
                            sample_order[sample][vcf_dictionary[vcf]]=i
                            i += 1
                elif "<ID=VARID," in line or "<ID=set," in line:
                    continue
                elif line[0] == line[1] and "=" in line:
                    if("ID=" in line and not "##contig=<ID=" in line):
                        field=line.split("=")[2].split(",")[0]
                        key= line.strip("#").split("=")[0]
                        if not key in header:
                            header[key]={}
                        header[key][field]=line
                    elif "##contig=<ID=" in line and not contigs:
                        header["CONTIGS"].append(line)
                    elif "##reference=" in line and not contigs:
                        reference=line
                    elif not "##source=" in line and not "##file" in line and not "##reference=" in line and not "##contig=<ID=" in line:
                        key=line.strip("#").split("=")[0]
                        if not key in subheader:
                                subheader[key]=line
            else:
                if header["CONTIGS"]:
                    contigs = True
                break

    #print the mandatory header lines in the correct order
    
    for entry in sorted(header["ALT"]):
        print(header["ALT"][entry].strip())
    del header["ALT"]
    for entry in sorted(header["INFO"]):
        print(header["INFO"][entry].strip())
    del header["INFO"]
    #print contigs according to the input order
    if reference != "":
        print reference.strip()
    for entry in header["CONTIGS"]:
        print(entry.strip())
        contigs_list.append(entry.strip())
    del header["CONTIGS"]
    for entry in sorted(header["FILTER"]):
        print(header["FILTER"][entry].strip())
    del header["FILTER"]

    all_formats=[]
    for entry in sorted(header["FORMAT"]):
        print(header["FORMAT"][entry].strip())
        

    del header["FORMAT"]

    #print the other lines in lexiographic order
    for key in sorted(header):
        for entry in sorted(header[key]):
            print(header[key][entry].strip())


    #print subheaders
    for entry in sorted(subheader):
        print(subheader[entry].strip())
    if not args.notag:
        print("##INFO=<ID=VARID,Number=1,Type=String,Description=\"The variant ID of merged samples\">")
        print("##INFO=<ID=set,Number=1,Type=String,Description=\"Source VCF for the merged record in SVDB\">")
    print("##svdbcmdline={}".format(" ".join(command_line)))
    sample_print_order={}

    if args.same_order:
        print first_vcf_header
    else:
        print( "\t".join(columns) )
    return (samples,sample_order,sample_print_order,contigs_list)

def main(args):
    variants={}
    #add all the variants into a dictionary
    i=0;
    vcf_list=args.vcf

    if vcf_list == 0:
        print("invalid input, either supply the vcf files, or add a number after each vcf name to asign the order manually")
        return(-1)
    else:
        vcf_dictionary={}
        priority_order = [] 
        if args.priority:
            priority_list=[]
            priority_dictionary={}
            vcf_dictionary={}
            for vcf in vcf_list:
                priority_dictionary[vcf.split(":")[-1]]=vcf.split(":")[0]
                vcf_dictionary[vcf.split(":")[0]] = vcf.split(":")[-1]
            for tag in args.priority.split(","):
                if tag in priority_dictionary:
                    priority_list.append(tag)

            if not len(priority_list) == len(priority_dictionary):
                print ("error tag/vcf mismatch, make sure that there is one tag per input vcf, or skip the --priority flag")
                return(-1)

            vcf_list=[]
            for tag in priority_list:
                vcf_list.append(priority_dictionary[tag])
            priority_order = priority_list

        for vcf in vcf_list:
            if not args.priority:
                vcf_dictionary[ vcf ]=vcf.split(".vcf")[0].split("/")[-1]
                priority_order.append(vcf.split(".vcf")[0].split("/")[-1])

            for line in open(vcf):
                if line[0] == "#":
                    pass
                else:
                    chrA,posA,chrB,posB,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
                    if not chrA in variants:
                        variants[chrA]=[]
                    if args.priority:
                        variants[chrA].append([chrB,event_type,posA,posB, vcf_dictionary[vcf],i,line.strip()])
                    else:
                        variants[chrA].append([chrB,event_type,posA,posB,vcf,i,line.strip()])
                    i+=1;

    samples,sample_order,sample_print_order,contigs=print_header(vcf_list,vcf_dictionary,args,sys.argv)
    to_be_printed=merge_vcf_module_cython.merge(variants,samples,sample_order,sample_print_order,priority_order,args)

    #use the contig order as defined in the header, or use lexiographic order
    if contigs:
        for i in range(0,len(contigs)):
            contig=contigs[i].split("##contig=<ID=")[-1].split(",length=")[0]
            if not contig in contigs:
                contigs[i]=contig

    else:
        contigs=sorted(to_be_printed.keys())

    #make sure all chromosomes were printed in the header
    for chromosome in to_be_printed:
        if not chromosome in contigs:
            contigs.append(chromosome)

    #print the variants in similar order as the input
    for chra in contigs:
        if not chra in to_be_printed:
            continue
        for variant in sorted(to_be_printed[chra],key = lambda x: int(x[1])):
            print("\t".join(variant).strip())
