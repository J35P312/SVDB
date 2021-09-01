from __future__ import absolute_import

from . import overlap_module


def retrieve_key(line, key):
    key += '='
    if key in line:
        item = line.strip().split(key)[-1].split(";")[0]
        if len(item) == len(line.strip()):
            return False
    return item


def determine_set_tag(priority_order, files):
    n_filtered = 0
    n_pass = 0

    filtered = []
    for sample in priority_order:
        file = files.get(sample, None)
        if file is None:
            continue
        if file.split('\t')[6] in ['PASS', '.']:
            n_pass += 1
        else:
            n_filtered += 1
    if n_pass == len(priority_order):
        return "Intersection"
    elif n_filtered == len(priority_order):
        return "FilteredInAll"
    else:
        for sample in priority_order:
            if sample not in files:
                continue
            elif files[sample].split('\t')[6] in ['PASS', '.']:
                filtered.append(sample)
            else:
                filtered.append("filterIn" + sample)
        return "-".join(filtered)


def merge_csq(info, csq):
    """Merge the csq fields of bnd variants"""
    var_csq = info.split("CSQ=")[-1].split(";")[0]
    csq.append(var_csq)
    effects = set([])
    for csq_field in csq:
        effects = effects | set(csq_field.split(","))
    CSQ = "CSQ=" + ",".join(list(effects))
    pre_CSQ = info.split("CSQ=")[0]
    post_CSQ = info.split("CSQ=")[-1].split(";")
    if len(post_CSQ) == 1:
        post_CSQ = ""
    else:
        post_CSQ = ";" + ";".join(post_CSQ[1:])

    return pre_CSQ + CSQ + post_CSQ


def sort_format_field(line, samples, sample_order, sample_print_order, priority_order, files, representing_file, args):
    format_columns = {}
    format_entries = []
    format_entry_length = []

    if not args.same_order:
        for input_file in priority_order:
            if input_file not in files:
                continue
            for sample in sorted(samples):
                if sample not in format_columns and input_file in sample_order[sample]:

                    vcf_line = files[input_file].strip().split("\t")
                    sample_position = sample_order[sample][input_file]
                    format_columns[sample] = {}
                    entries = vcf_line[8].split(":")
                    sample_entries = vcf_line[9 + sample_position].split(":")

                    for i, entry in enumerate(entries):
                        format_columns[sample][entry] = sample_entries[i]
                        if entry not in format_entries:
                            n = sample_entries[i].count(",")
                            format_entries.append(entry)
                            format_entry_length.append(n)

        if len(line) > 8:
            format_string = []
            line[8] = ":".join(format_entries)
            del line[9:]
            for sample in samples:
                format_string = []
                for entry in format_entries:
                    j = 0
                    if sample in format_columns:
                        if entry in format_columns[sample]:
                            format_string.append(format_columns[sample][entry])
                        elif entry == "GT":
                            format_string.append("./.")
                        else:
                            sub_entry = []
                            for i in range(format_entry_length[j] + 1):
                                sub_entry.append(".")
                            format_string.append(",".join(sub_entry))
                    else:
                        if entry == "GT":
                            format_string.append("./.")
                        else:
                            sub_entry = []
                            for i in range(format_entry_length[j] + 1):
                                sub_entry.append(".")
                            format_string.append(",".join(sub_entry))
                    j += 1

                line.append(":".join(format_string))

    # generate a union of the info fields
    info_union = []
    tags_in_info = []
    for input_file in priority_order:
        if input_file not in files:
            continue
        INFO = files[input_file].strip().split("\t")[7]
        INFO_content = INFO.split(";")

        for content in INFO_content:
            tag = content.split("=")[0]
            if tag not in tags_in_info:
                tags_in_info.append(tag)
                info_union.append(content)

    new_info = ";".join(info_union)
    line[7] = new_info
    return line


def merge(variants, samples, sample_order, sample_print_order, priority_order, args):
    overlap_param = args.overlap
    bnd_distance = args.bnd_distance
    ins_distance=args.ins_distance
    no_intra = args.no_intra
    no_var = args.no_var
    pass_only = args.pass_only

    # search for similar variants
    to_be_printed = {}
    for chrA in variants:
        analysed_variants = set([])
        for i in range(len(variants[chrA])):
            if i in analysed_variants:
                continue

            merge = []
            csq = []

            files = {}
            for j in range(i + 1, len(variants[chrA])):
                if j in analysed_variants:
                    continue

                # if the pass_only option is chosen, only variants marked PASS will be merged
                if pass_only:
                    filter_tag = variants[chrA][i][-1].split("\t")[6]
                    if filter_tag not in ['PASS', '.']:
                        break

                # only treat varints on the same pair of chromosomes
                if not variants[chrA][i][0] == variants[chrA][j][0]:
                    continue

                # if the pass_only option is chosen, only variants marked PASS will be merged
                if pass_only:
                    filter_tag = variants[chrA][j][-1].split("\t")[6]
                    if filter_tag not in ['PASS', '.']:
                        continue

                # dont merge variants of different type
                if variants[chrA][i][1] != variants[chrA][j][1] and not no_var:
                    continue

                # if no_intra is chosen, variants may only be merged if they belong to different input files
                if no_intra and variants[chrA][i][-3] == variants[chrA][j][-3]:
                    continue

                if variants[chrA][i][0] == "INS":
                    overlap, match = overlap_module.variant_overlap(
                        chrA, variants[chrA][i][0], variants[chrA][i][2], variants[chrA][i][3], variants[chrA][j][2], variants[chrA][j][3], -1, ins_distance)

                else:
                    overlap, match = overlap_module.variant_overlap(
                        chrA, variants[chrA][i][0], variants[chrA][i][2], variants[chrA][i][3], variants[chrA][j][2], variants[chrA][j][3], overlap_param, bnd_distance)


                if match:
                    # add similar variants to the merge list and remove them
                    if args.priority:
                        files[variants[chrA][j][-3]] = variants[chrA][j][-1]
                        merge.append(
                            variants[chrA][j][-1].split("\t")[2].replace(";", "_") + ":" + variants[chrA][j][-3])
                    else:
                        files[variants[chrA][j]
                              [-3].split(".vcf")[0].split("/")[-1]] = variants[chrA][j][-1]
                        merge.append(variants[chrA][j][-1].split("\t")[2].replace(
                            ";", "_") + ":" + variants[chrA][j][-3].split(".vcf")[0].split("/")[-1])

                    if variants[chrA][i][0] != chrA and "CSQ=" in variants[chrA][j][-1]:
                        info = variants[chrA][j][-1].split("\t")[7]
                        csq.append(info.split("CSQ=")[-1].split(";")[0])
                    analysed_variants.add(j)

            line = variants[chrA][i][-1].split("\t")

            if csq:
                line[7] = merge_csq(line[7], csq)
            if not line[0] in to_be_printed:
                to_be_printed[line[0]] = []

            if args.priority:
                files[variants[chrA][i][-3]] = "\t".join(line)
                representing_file = variants[chrA][i][-3]
            else:
                files[variants[chrA][i]
                      [-3].split(".vcf")[0].split("/")[-1]] = "\t".join(line)
                representing_file = variants[chrA][i][-3].split(".vcf")[
                    0].split("/")[-1]
            line = sort_format_field(
                line, samples, sample_order, sample_print_order, priority_order, files, representing_file, args)
            if merge and not args.notag:
                line[7] += ";VARID=" + "|".join(merge)
                line[2] += ":{}|".format(variants[chrA][i][-3].split(".vcf")
                                         [0].split("/")[-1]) + "|".join(merge)
            if not args.notag:
                set_tag = determine_set_tag(priority_order, files)
                line[7] += ";set={}".format(set_tag)
            to_be_printed[line[0]].append(line)

            analysed_variants.add(i)

    return to_be_printed
