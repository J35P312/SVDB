from __future__ import absolute_import

import glob
import gzip
import os
import sqlite3

from . import readVCF


def populate_db(args):
    sample_IDs = []
    conn = sqlite3.connect(args.db + ".db")
    c = conn.cursor()
    tableListQuery = "SELECT name FROM sqlite_master WHERE type=\'table\'"
    c.execute(tableListQuery)
    tables = map(lambda t: t[0], c.fetchall())

    idx = 0
    if "SVDB" not in tables:
        A = "CREATE TABLE SVDB (var TEXT,chrA TEXT, chrB TEXT,posA INT,ci_A_lower INT,ci_A_upper INT,posB INT,ci_B_lower INT,ci_B_upper INT, sample TEXT, idx INT)"
        c.execute(A)

    else:
        try:
            A = "DROP INDEX SV "
            c.execute(A)
        except Exception:
            pass
        try:
            A = "DROP INDEX IDX"
            c.execute(A)
        except Exception:
            pass
        try:
            A = "DROP INDEX CHR"
            c.execute(A)
        except Exception:
            pass

        A = 'SELECT DISTINCT sample FROM SVDB'
        for sample in c.execute(A):
            sample_IDs.append(sample)

        if sample_IDs:
            A = "SELECT MAX(idx) FROM SVDB"
            number_of_variants = 0
            for hit in c.execute(A):
                number_of_variants = int(hit[0])
            idx = 1 + number_of_variants

    # populate the tables
    for vcf in args.files:
        sample_name = vcf.split("/")[-1].split(".vcf")[0]
        sample_name = sample_name.replace(".", "_")
        sample_IDs.append(sample_name)
        A = 'SELECT sample FROM SVDB WHERE sample == \'{}\' '.format(
            sample_name)
        hits = []
        for hit in c.execute(A):
            hits.append(hit)

        if hits:
            continue

        if not os.path.exists(vcf):
            print("error: unnable to open {}".format(vcf))
            continue

        var = []
        sample_names = []

        if vcf.endswith(".gz"):
            f = gzip.open(vcf, "rt")
        else:
            f = open(vcf, "rt")

        for line in f:
            if line.startswith("#"):
                if "CHROM" in line:
                    content = line.strip().split()
                    if len(content) > 9:
                        sample_names = content[9:]
                continue

            if not len(line.strip()):
                continue

            chrA, posA, chrB, posB, event_type, INFO, FORMAT = readVCF.readVCFLine(
                line)
            if args.passonly:
                FILTER = line.split("\t")[6]
                if not (FILTER == "PASS" or FILTER == "."):
                    continue

            ci_A_lower = 0
            ci_A_upper = 0
            ci_B_lower = 0
            ci_B_upper = 0
            if "CIPOS" in INFO:
                ci = INFO["CIPOS"].split(",")
                if len(ci) > 1:
                    ci_A_lower = abs(int(ci[0]))
                    ci_A_upper = abs(int(ci[1]))
                    ci_B_lower = abs(int(ci[0]))
                    ci_B_upper = abs(int(ci[1]))
                else:
                    ci_A_lower = abs(int(ci[0]))
                    ci_A_upper = abs(int(ci[0]))
                    ci_B_lower = abs(int(ci[0]))
                    ci_B_upper = abs(int(ci[0]))

            if "CIEND" in INFO:
                ci = INFO["CIEND"].split(",")
                if len(ci) > 1:
                    ci_B_lower = abs(int(ci[0]))
                    ci_B_upper = abs(int(ci[1]))
                else:
                    ci_B_lower = abs(int(ci[0]))
                    ci_B_upper = abs(int(ci[0]))

            if "GT" not in FORMAT or not len(sample_names):
                var.append((event_type, chrA, chrB, posA, ci_A_lower,
                            ci_A_upper, posB, ci_B_lower, ci_B_upper, sample_name, idx))
                idx += 1
            else:
                sample_index = 0
                for genotype in FORMAT["GT"]:
                    if genotype not in ["0/0", "./."]:
                        var.append((event_type, chrA, chrB, posA, ci_A_lower, ci_A_upper,
                                    posB, ci_B_lower, ci_B_upper, sample_names[sample_index], idx))
                        idx += 1
                    sample_index += 1

            # insert EVERYTHING into the database, the user may then query it in different ways(at least until the DB gets to large to function properly)
        if var:
            c.executemany(
                'INSERT INTO SVDB VALUES (?,?,?,?,?,?,?,?,?,?,?)', var)

    A = "CREATE INDEX SV ON SVDB (var,chrA, chrB, posA,posA,posB,posB)"
    c.execute(A)
    A = "CREATE INDEX IDX ON SVDB (idx)"
    c.execute(A)
    A = "CREATE INDEX CHR ON SVDB (chrA, chrB)"
    c.execute(A)

    conn.commit()
    conn.close()
    return(sample_IDs)


def main(args):
    args.db = args.prefix
    if not args.files and args.folder:
        args.files = glob.glob(os.path.join(args.folder, "*.vcf")) + glob.glob(os.path.join(args.folder, "*.vcf.gz"))
    populate_db(args)
