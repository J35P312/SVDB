import readVCF
import SVDB_merge_vcf_module_cython

def check_input(vcf_list):
    number_test =[j for j in vcf_list if j.isdigit()]
    #there are numbers in the list
    if number_test:
        if len(number_test) == len(vcf_list)/2:
            number_list=[int(j) for j in vcf_list[1::2] if j.isdigit()]

            expected_numbers=range(1,len(vcf_list)/2+1)
            if len(expected_numbers) == len(number_list) and set(number_list) == set(expected_numbers):
                vcf_dictionary=dict(zip(vcf_list[1::2], vcf_list[0::2]))
                vcf_list=[]
                for order in sorted(vcf_dictionary):
                    vcf_list.append(vcf_dictionary[order])
                return(vcf_list)
        return(0);
    else:
        return(vcf_list);

def print_header(vcf_list):

    header={};
    header["ALT"]={}
    header["INFO"]={}
    header["FILTER"]={}
    header["FORMAT"]={}

    subheader={}
    columns=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    
    print("##fileformat=VCFv4.1")
    print("##source=MergeVCF")
    for vcf in vcf_list:
        with open(vcf) as text_file:
            for line in text_file:
                if(line[0] == "#"):
                    if("#CHROM\tPOS" in line):
                        vcf_columns=line.strip().split("\t")
                        for column in vcf_columns:
                            if not column in columns and len(vcf_columns) > len(columns):
                                columns.append(column.strip())
                        break
                    
                    elif line[0] == line[1] and line[0] == "#" and "=" in line:
                        if("ID=" in line):
                            field=line.split("=")[2].split(",")[0]
                            key= line.strip("#").split("=")[0]
                            if not key in header:
                                header[key]={field:line}
                            else:
                                header[key].update({field:line})
                        elif not "##source=" in line and not "##file" in line and not "##reference=" in line:
                            key= line.strip("#").split("=")[0]
                            if not key in subheader:
                                subheader.update({key:line})
                else:
                    break

    #print the mandatory header lines in the correct order
    for entry in sorted(header["ALT"]):
        print(header["ALT"][entry].strip())
    del header["ALT"]
    for entry in sorted(header["INFO"]):
        print(header["INFO"][entry].strip())
    del header["INFO"]
    for entry in sorted(header["FILTER"]):
        print(header["FILTER"][entry].strip())
    del header["FILTER"]

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
    print("##INFO=<ID=VARID,Number=1,Type=String,Description=\"The variant ID of merged samples\">")
    print( "\t".join(columns) )
    return("nada")

def main(args):
    variants={}
    #add all the variants into a dictionary
    i=0;
    vcf_list=check_input(args.vcf)
    if vcf_list == 0:
        print("invalid input, either supply the vcf files, or add a number after each vcf name to asign the order manually")
        return(-1)
    else:
        for vcf in vcf_list:
            for line in open(vcf):
                if line[0] == "#":
                    pass
                else:
                    chrA,posA,chrB,posB,event_type,INFO,FORMAT =readVCF.readVCFLine(line);
                    if not chrA in variants:
                        variants[chrA]=[]
                    variants[chrA].append([chrB,event_type,posA,posB,vcf,i,line.strip()])
                    i+=1;

    print_header(vcf_list)

    to_be_printed=SVDB_merge_vcf_module_cython.merge(variants,args.ci,args.overlap,args.bnd_distance,args.no_intra)
    #print the variants in similar order as the input
    for chra in sorted(to_be_printed):
        for variant in sorted(to_be_printed[chra],key = lambda x: int(x[1])):
            print("\t".join(variant))
