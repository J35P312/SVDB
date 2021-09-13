import re


# TODO: Should be part of a VCF class
def readVCFLine(line):
    if line[0].startswith("#"):
        return None

    variation = line.strip().split("\t")
    event_type = ""
    chrA = variation[0].replace("chr", "").replace("Chr", "").replace("CHR", "")
    posA = int(variation[1])
    posB = 0

    description = {}
    INFO = variation[7].split(";")
    for tag in INFO:
        tag = tag.split("=")
        if(len(tag) > 1):
            description[tag[0]] = tag[1]

    format = {}
    format_keys = {}
    if len(variation) > 8:
        format_string = variation[8].split(":")

        i = 0
        for key in format_string:
            format_keys[i] = key
            format[key] = []
            i += 1
        format_fields = variation[9:]
        for sample in format_fields:
            i = 0
            format_string = sample.split(":")
            for i in range(0, len(format_keys)):
                format[format_keys[i]].append(format_string[i])
                i += 1

    # Delly translocations
    if "TRA" in variation[4]:
        event_type = "BND"
        chrB = description["CHR2"]
        posB = int(description["END"])
        if chrA > chrB:
            chrT = chrA
            chrA = chrB
            chrB = chrT
            posB = posA

    # intrachromosomal variant
    elif "]" not in variation[4] and "[" not in variation[4]:

        chrB = chrA
        if "END" in description:
            posB = int(description["END"])

        elif "SVLEN" in description:
            posB = posA + abs(int(description["SVLEN"]))
        else:
            posB = posA
        # sometimes the fermikit intra chromosomal events are inverted i.e the end pos is a lower position than the start pos
        if posB < posA:
            tmp = posB
            posB = posA
            posA = tmp
        event_type = variation[4].strip("<").rstrip(">")

        if "<" in variation[4] and ">" in variation[4]:
            if "DUP" in event_type:
                event_type = "DUP"
        else:
            if "SVTYPE" in description:
                event_type = description["SVTYPE"]
        #treat the insertion as single points
        if "INS" in event_type:
            posA=int(variation[1])
            posB=int(variation[1])



    # if the variant is given as a breakpoint, it is stored as a precise variant in the db
    else:
        B = variation[4]
        combinations={"[]":"]", "[[":"[", "]]":"]", "][":"["}
        
        for c in combinations:
            if B.startswith(c):
                B=B.replace(c,combinations[c])

        B = re.split("[],[]", B)
        chr_and_pos = B[1]
        chrB = ":".join(chr_and_pos.split(":")[:-1]).replace("chr", "").replace("Chr", "").replace("CHR", "")
        posB = int(chr_and_pos.split(":")[-1])
        if chrA > chrB:
            chrT = chrA
            chrA = chrB
            chrB = chrT
        posA, posB = posB, posA
        event_type = "BND"

    return chrA, posA, chrB, posB, event_type, description, format
