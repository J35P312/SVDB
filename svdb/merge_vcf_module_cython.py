from . import overlap_module


def sanitize_id(s: str) -> str:
    """Replace VCF-unsafe delimiters in a variant ID with underscores."""
    return s.replace(";", "_").replace(":", "_").replace("|", "_")


def format_tag(var_id: str, value: str) -> str:
    """Format a per-variant tag entry as 'sanitized_id|value'."""
    return f"{sanitize_id(var_id)}|{value}"


def retrieve_key(line, key):
    key += '='
    if key not in line:
        return False
    if f";{key}" in line:
        return line.strip().split(f";{key}")[-1].split(";")[0].split("\t")[0]
    if f"\t{key}" in line:
        return line.strip().split(f"\t{key}")[-1].split(";")[0].split("\t")[0]
    return False

#Check if no merging should occur
def skip_variant(chrA,chrB,type_A,type_B,vcf_line_A,vcf_line_B,pass_only,current_variant,analysed_variants,no_var):
    #The variant is already clustered/analysed
    if current_variant in analysed_variants:
       return True

    # only treat variants on the same pair of chromosomes
    if chrA != chrB:
       return True

    # dont merge variants of different type
    if type_A != type_B and not no_var:
       return True

    # if the pass_only option is chosen, only variants marked PASS will be merged
    if pass_only:
       filter_tag = vcf_line_B[6]
       if filter_tag not in ['PASS', '.']:
           return True

    return False


#Collect SAMPLE columns from all merged variants:
def collect_sample(vcf_line,samples,sample_order,f):
    variant=vcf_line[2].replace(";","_").replace(":","_").replace("|","_")
    sample_data=[variant]

    for sample in samples:
        if sample not in sample_order:
           continue
        if f not in sample_order[sample]:
           continue
        sample_position = sample_order[sample][f]

        entries = vcf_line[8].split(":")
        sample_entries = vcf_line[9 + sample_position].split(":")
        sample_data.append(sample)
        for i, entry in enumerate(entries):
            sample_data.append(f"{entry}:{sample_entries[i]}" )

    return "|".join(sample_data).replace(",",":")

#collect INFO from all merged_variants
def collect_info(vcf_line):
    INFO = vcf_line[7]
    INFO_content = INFO.split(";")
    variant=vcf_line[2].replace(";","_").replace(":","_").replace("|","_")
    all_info=[variant]

    for content in INFO_content:
        if ":" not in content and "|" not in content:
           all_info.append( content.replace("=",":").replace(",",":") )

    return "|".join(all_info)

#create a GATK-like set of all merged variants
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

#merge csq field of merged variants (for instance when merging BNDs)
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


def sort_format_field(line, samples, sample_order, priority_order, files, args):
    format_columns = {}
    format_entries = []
    format_entry_length = []

    if not args.same_order:
        sorted_samples = sorted(samples)  # constant across all input_files — compute once
        for input_file in priority_order:
            if input_file not in files:
                continue
            for sample in sorted_samples:
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
    tags_in_info = set()  # set for O(1) membership test vs O(n) list

    # tags only to be copied from the file with highest priority (to avoid problems in downstream analyses
    blacklist=set(["SVLEN","END","SVTYPE"])

    first=True
    for input_file in priority_order:
        if input_file not in files:
            continue

        INFO = files[input_file].strip().split("\t")[7]
        INFO_content = INFO.split(";")

        for content in INFO_content:
            tag = content.split("=")[0]

            if not first and tag in blacklist:
                continue

            if tag not in tags_in_info:
                tags_in_info.add(tag)
                info_union.append(content)

        first=False

    new_info = ";".join(info_union)
    line[7] = new_info
    return line


def merge(variants, samples, sample_order, priority_order, args):
    overlap_param = args.overlap
    bnd_distance = args.bnd_distance
    ins_distance=args.ins_distance
    no_intra = args.no_intra
    no_var = args.no_var
    pass_only = args.pass_only

    # search for similar variants
    to_be_printed = {}
    for chrA in variants:
        # Sort by posA so the positional early-exit break below is valid.
        # Input variants[chrA] is built by concatenating multiple VCF files
        # which may themselves be unsorted, so we normalise here.
        # O(n log n) is negligible compared to the O(n²) loop it enables.
        variants[chrA].sort(key=lambda v: v.posA)
        analysed_variants = set([])
        for i in range(len(variants[chrA])):
            if i in analysed_variants:
                continue

            merge = []
            #keep track of all csq of merged variants
            csq = []

            #Keep track of all files in the cluster
            files = {}
            #keep track of the FILTER column of all merged files
            filters_tag = {}
            #keep track of SAMPLEs columns of all merged files
            samples_tag = {}
            #keep track of INFO column of all merged files
            info_tag = {}
            #quality
            qual_tag = {}
            #pos
            pos_tag = {}
            #chrom
            chrom_tag ={}

            # Cache variant_i attributes — avoids repeated attribute lookups inside the O(n²) j-loop
            var_i = variants[chrA][i]
            source_i = var_i.source
            chrB_i = var_i.chrB
            type_i = var_i.event_type
            posA_i = var_i.posA
            posB_i = var_i.posB
            insertion_i = var_i.is_insertion()  # pre-compute once; was called per-j-iteration

            if args.priority:
               id = source_i
            else:
               id = source_i.split(".vcf")[0].split("/")[-1]

            vcf_line_A=var_i.raw_line.strip().split("\t")
            chrom_tag[id]=[ format_tag(vcf_line_A[2], vcf_line_A[0]) ]
            pos_tag[id]=[ format_tag(vcf_line_A[2], vcf_line_A[1]) ]
            qual_tag[id]=[ format_tag(vcf_line_A[2], vcf_line_A[5]) ]
            filters_tag[id]=[ format_tag(vcf_line_A[2], vcf_line_A[6]) ]
            samples_tag[id]=[collect_sample( vcf_line_A ,samples,sample_order,id)]
            info_tag[id]=[collect_info(vcf_line_A)]

            # pass_only: variant A must pass filter to merge others into it.
            # This condition is constant across all j — evaluate once outside the loop.
            a_can_merge = (not pass_only) or vcf_line_A[6] in ('PASS', '.')

            for j in range(i + 1, len(variants[chrA])) if a_can_merge else []:
                # --- Cheap early-exit checks (no string split needed) ---
                if j in analysed_variants:
                    continue

                var_j = variants[chrA][j]

                # VCF files are sorted by (chrom, posA).  Once var_j.posA is
                # further than bnd_distance ahead of posA_i, no subsequent j
                # can overlap either — break out of the inner loop entirely.
                # BND/interchromosomal pairs skip this check (chrA != chrB_i).
                if chrA == chrB_i and var_j.posA - posA_i > bnd_distance:
                    break

                if chrB_i != var_j.chrB:
                    continue

                if type_i != var_j.event_type and not no_var:
                    continue

                # pass_only: need to inspect the filter field — split early only when required
                if pass_only:
                    vcf_line_B = var_j.raw_line.strip().split("\t")
                    if vcf_line_B[6] not in ('PASS', '.'):
                        continue

                # if no_intra is chosen, variants may only be merged if they belong to different input files
                if no_intra and source_i == var_j.source:
                    continue

                if insertion_i:
                    overlap, match = overlap_module.variant_overlap(
                        chrA, chrB_i, posA_i, posB_i, var_j.posA, var_j.posB, -1, ins_distance)
                else:
                    overlap, match = overlap_module.variant_overlap(
                        chrA, chrB_i, posA_i, posB_i, var_j.posA, var_j.posB, overlap_param, bnd_distance)

                if match:
                    # Split only on confirmed match (pass_only case already split above)
                    if not pass_only:
                        vcf_line_B = var_j.raw_line.strip().split("\t")

                    # add similar variants to the merge list and remove them
                    if args.priority:
                        match_id = var_j.source
                    else:
                        match_id = var_j.source.split(".vcf")[0].split("/")[-1]

                    files[match_id] = var_j.raw_line
                    merge.append(sanitize_id(vcf_line_B[2]) + ":" + match_id)
                    if match_id not in filters_tag:
                        filters_tag[match_id]=[]
                        samples_tag[match_id]=[]
                        info_tag[match_id]=[]
                        chrom_tag[match_id]=[]
                        pos_tag[match_id]=[]
                        qual_tag[match_id]=[]

                    chrom_tag[match_id].append(format_tag(vcf_line_B[2], vcf_line_B[0]))
                    pos_tag[match_id].append(format_tag(vcf_line_B[2], vcf_line_B[1]))
                    qual_tag[match_id].append(format_tag(vcf_line_B[2], vcf_line_B[5]))
                    filters_tag[match_id].append(format_tag(vcf_line_B[2], vcf_line_B[6]))

                    samples_tag[match_id].append(collect_sample(vcf_line_B, samples, sample_order, match_id))
                    info_tag[match_id].append(collect_info(vcf_line_B))

                    if chrB_i != chrA and "CSQ=" in var_j.raw_line:
                        info = vcf_line_B[7]
                        csq.append(info.split("CSQ=")[-1].split(";")[0])
                    analysed_variants.add(j)

            line = vcf_line_A

            if csq:
                line[7] = merge_csq(line[7], csq)
            if line[0] not in to_be_printed:
                to_be_printed[line[0]] = []

            if args.priority:
                files[source_i] = "\t".join(line)
            else:
                files[source_i.split(".vcf")[0].split("/")[-1]] = "\t".join(line)
            line = sort_format_field(
                line, samples, sample_order, priority_order, files, args)
            if merge and not args.notag:
                line[7] += ";VARID=" + "|".join(merge)
                line[2] += ":{}|".format(source_i.split(".vcf")[0].split("/")[-1]) + "|".join(merge)
            if not args.notag:
                set_tag = determine_set_tag(priority_order, files)
                line[7] += f";set={set_tag}"
                line[7] += f";FOUNDBY={len(set(filters_tag.keys()))}"

            #add chrom information of all merged variants
            callers=[]
            for tag in filters_tag:
                line[7]+=";{}_CHROM={}".format(tag,",".join(chrom_tag[tag]))
                callers.append(tag)
            #add pos information of all merged variants
            for tag in filters_tag:
                line[7]+=";{}_POS={}".format(tag,",".join(pos_tag[tag]))
            #add qual information of all merged variants
            for tag in filters_tag:
                line[7]+=";{}_QUAL={}".format(tag,",".join(qual_tag[tag]))
            #add filter of all merged variants
            for tag in filters_tag:
                line[7]+=";{}_FILTERS={}".format(tag,",".join(filters_tag[tag]).replace(";",","))
            #add samples information for all merged variants
            for tag in samples_tag:
                line[7]+=";{}_SAMPLE={}".format(tag,",".join(samples_tag[tag]))
            #add info column for all merged variants
            for tag in samples_tag:
                line[7]+=";{}_INFO={}".format(tag,",".join(info_tag[tag]))

            if not args.notag:
                line[7]+=";svdb_origin={}".format("|".join(callers))

            sup_vec=[]
            for c in sorted(priority_order):
                if c in filters_tag:
                    sup_vec.append("1")
                else:
                    sup_vec.append("0")

            line[7]+=";SUPP_VEC={}".format("".join(sup_vec))

            to_be_printed[line[0]].append(line)

            analysed_variants.add(i)

    return to_be_printed
