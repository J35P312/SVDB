import logging
import sys

from . import merge_vcf_module_cython, read_vcf, vcf_utils
from .models import MergeVariant

logger = logging.getLogger(__name__)


def build_header(vcf_list, vcf_dictionary, args, command_line):
    """Build the VCF merge header string and collect sample metadata.

    Returns (header_string, samples, sample_order, contigs_list).
    No I/O side effects — the caller is responsible for writing the header.
    """
    header = {"ALT": {},
              "INFO": {},
              "FILTER": {},
              "FORMAT": {},
              "CONTIGS": []}
    contigs = False
    reference = ""
    subheader = {}
    contigs_list = []
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    sample_order = {}
    lines_out = []
    lines_out.append("##fileformat=VCFv4.1")
    lines_out.append("##source=MergeVCF")
    lines_out.append("##SVDB_version={} cmd=\"{}\"".format(args.version, " ".join(sys.argv)))
    samples = []
    first = True
    first_vcf_header = ""
    for vcf in vcf_list:
        with vcf_utils.open_vcf(vcf) as lines:
            for line in lines:
                if line.startswith('#'):
                    if "#CHROM\tPOS" in line:
                        if first:
                            first_vcf_header = line.strip()
                            first = False
                        vcf_columns = line.strip().split("\t")
                        for column in vcf_columns:
                            if column not in columns:
                                columns.append(column)
                                if column not in sample_order and column != "FORMAT":
                                    sample_order[column] = {}
                                    samples.append(column)

                        if len(vcf_columns) > 8 and not args.same_order:
                            for i, sample in enumerate(vcf_columns[9:]):
                                sample_order[sample][vcf_dictionary[vcf]] = i
                    elif "<ID=VARID," in line or "<ID=set," in line:
                        continue
                    elif line[0] == line[1] and "=" in line:
                        if("ID=" in line and "##contig=<ID=" not in line):
                            field = line.split("=")[2].split(",")[0]
                            key = line.strip("#").split("=")[0]
                            if key not in header:
                                header[key] = {}
                            header[key][field] = line
                        elif "##contig=<ID=" in line and not contigs:
                            header["CONTIGS"].append(line)
                        elif "##reference=" in line and not contigs:
                            reference = line
                        elif "##source=" not in line and "##file" not in line and "##reference=" not in line and "##contig=<ID=" not in line:
                            key = line.strip("#").split("=")[0]
                            if key not in subheader:
                                subheader[key] = line
                else:
                    if header["CONTIGS"]:
                        contigs = True
                    break

    # mandatory header lines in the correct order
    for entry in sorted(header["ALT"]):
        lines_out.append(header["ALT"][entry].strip())
    del header["ALT"]
    for entry in sorted(header["INFO"]):
        lines_out.append(header["INFO"][entry].strip())
    del header["INFO"]

    supp_vector_order = []
    for vcf in vcf_dictionary:
        lines_out.append("##INFO=<ID={}_INFO,Number=.,Type=String,Description=\"pipe separated list of all details in the INFO column of file {}\">".format(vcf_dictionary[vcf], vcf_dictionary[vcf]))
        lines_out.append("##INFO=<ID={}_SAMPLE,Number=.,Type=String,Description=\"pipe separated list of all details in the SAMPLEs column of file {}\">".format(vcf_dictionary[vcf], vcf_dictionary[vcf]))
        lines_out.append("##INFO=<ID={}_CHROM,Number=.,Type=String,Description=\"pipe separated list of all details in the CHROM column of file {}\">".format(vcf_dictionary[vcf], vcf_dictionary[vcf]))
        lines_out.append("##INFO=<ID={}_POS,Number=.,Type=String,Description=\"pipe separated list of all details in the POS column of file {}\">".format(vcf_dictionary[vcf], vcf_dictionary[vcf]))
        lines_out.append("##INFO=<ID={}_QUAL,Number=.,Type=String,Description=\"pipe separated list of all details in the QUAL column of file {}\">".format(vcf_dictionary[vcf], vcf_dictionary[vcf]))
        lines_out.append("##INFO=<ID={}_FILTERS,Number=.,Type=String,Description=\"pipe separated list of all details in the FILTER column of file {}\">".format(vcf_dictionary[vcf], vcf_dictionary[vcf]))
        supp_vector_order.append(vcf_dictionary[vcf])

    lines_out.append("##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description=\"Vector of supporting callers/files (order: {}).\">".format(" ".join(supp_vector_order)))

    # contigs in input order
    if reference != "":
        lines_out.append(reference.strip())
    for entry in header["CONTIGS"]:
        lines_out.append(entry.strip())
        contigs_list.append(entry.strip())
    del header["CONTIGS"]
    for entry in sorted(header["FILTER"]):
        lines_out.append(header["FILTER"][entry].strip())
    del header["FILTER"]

    for entry in sorted(header["FORMAT"]):
        lines_out.append(header["FORMAT"][entry].strip())
    del header["FORMAT"]

    # other lines in lexicographic order
    for key in sorted(header):
        for entry in sorted(header[key]):
            lines_out.append(header[key][entry].strip())

    # subheaders
    for entry in sorted(subheader):
        lines_out.append(subheader[entry].strip())
    if not args.notag:
        lines_out.append("##INFO=<ID=VARID,Number=1,Type=String,Description=\"The variant ID of merged samples\">")
        lines_out.append("##INFO=<ID=FOUNDBY,Number=1,Type=Integer,Description=\"The number of files containing the variant\">")
        lines_out.append("##INFO=<ID=set,Number=1,Type=String,Description=\"Source VCF for the merged record in SVDB\">")
        lines_out.append("##INFO=<ID=svdb_origin,Number=1,Type=String,Description=\"pipe separated list of the VCF for the merged record in SVDB\">")

    lines_out.append("##svdbcmdline={}".format(" ".join(command_line)))

    if args.same_order:
        lines_out.append(first_vcf_header)
    else:
        lines_out.append("\t".join(columns))

    header_string = "\n".join(lines_out)
    return header_string, samples, sample_order, contigs_list


def print_header(vcf_list, vcf_dictionary, args, command_line):
    """Build and print the VCF merge header; return sample metadata."""
    header_string, samples, sample_order, contigs_list = build_header(
        vcf_list, vcf_dictionary, args, command_line
    )
    print(header_string)
    return samples, sample_order, contigs_list


def main(args):
    variants = {}
    # add all the variants into a dictionary
    i = 0
    vcf_list = args.vcf

    if not vcf_list:
        logger.error("invalid input: supply vcf files, or add a number after each vcf name to assign the order manually")
        return -1
    else:
        vcf_dictionary = {}
        priority_order = []
        if args.priority:
            priority_list = []
            priority_dictionary = {}
            vcf_dictionary = {}
            for vcf in vcf_list:
                priority_dictionary[vcf.split(":")[-1]] = vcf.split(":")[0]
                vcf_dictionary[vcf.split(":")[0]] = vcf.split(":")[-1]
            for tag in args.priority.split(","):
                if tag in priority_dictionary:
                    priority_list.append(tag)

            if len(priority_list) != len(priority_dictionary):
                logger.error("tag/vcf mismatch: ensure one tag per input vcf, or omit --priority")
                return -1

            vcf_list = []
            for tag in priority_list:
                vcf_list.append(priority_dictionary[tag])
            priority_order = priority_list

        for vcf in vcf_list:
            if not args.priority:
                vcf_dictionary[vcf] = vcf.split(".vcf")[0].split("/")[-1]
                priority_order.append(vcf.split(".vcf")[0].split("/")[-1])

            with vcf_utils.open_vcf(vcf) as lines:
                for line in lines:
                    if line.startswith('#'):
                        continue
                    else:
                        v = read_vcf.readVCFLine(line)
                        if v.chrA not in variants:
                            variants[v.chrA] = []
                        if args.priority:
                            variants[v.chrA].append(MergeVariant(v.chrB, v.event_type, v.posA, v.posB, vcf_dictionary[vcf], i, line.strip()))
                        else:
                            variants[v.chrA].append(MergeVariant(v.chrB, v.event_type, v.posA, v.posB, vcf, i, line.strip()))
                        i += 1

    samples, sample_order, contigs = print_header(vcf_list, vcf_dictionary, args, sys.argv)
    to_be_printed = merge_vcf_module_cython.merge(variants, samples, sample_order, priority_order, args)

    # use the contig order as defined in the header, or use lexiographic order
    if contigs:
        for i, contig in enumerate(contigs):
            contig = contig.split("##contig=<ID=")[-1].split(",length=")[0]
            if contig not in contigs:
                contigs[i] = contig
    else:
        contigs = sorted(to_be_printed.keys())

    # make sure all chromosomes were printed in the header
    for chromosome in to_be_printed:
        if chromosome not in contigs:
            contigs.append(chromosome)

    # print the variants in similar order as the input
    for chra in contigs:
        if chra not in to_be_printed:
            continue
        for variant in sorted(to_be_printed[chra], key=lambda x: int(x[1])):
            print("\t".join(variant).strip())
