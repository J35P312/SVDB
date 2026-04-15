import glob
import gzip
import os
from pathlib import Path

from . import database, readVCF


def populate_db(args):
    db = database.DB(args.db)
    tables = db.tables

    idx = 0
    if "SVDB" not in tables:
        db.create(database.CREATE_TABLE_SQL)
        sample_IDs = []
    else:
        db.drop("DROP INDEX SV")
        db.drop("DROP INDEX IDX")
        db.drop("DROP INDEX CHR")

        sample_IDs = db.sample_ids
        if sample_IDs:
            idx = 1 + int(db.query("SELECT MAX(idx) FROM SVDB")[0][0])

    # populate the tables
    for vcf in args.files:
        sample_name = Path(vcf).stem.replace(".", "_")
        sample_IDs.append(sample_name)
        A = 'SELECT sample FROM SVDB WHERE sample == \'{}\' '.format(sample_name)
        hits = [hit for hit in db.query(A)]
        if hits:
            continue
        if not os.path.exists(vcf):
            print("error: unnable to open {}".format(vcf))
            continue

        var = []
        sample_names = []

        # TODO: Move this into a VCF class
        opener = gzip.open if vcf.endswith('.vcf.gz') else open
        with opener(vcf, 'rt') as lines:
            for line in lines:
                if line.startswith("#"):
                    if "CHROM" in line:
                        content = line.strip().split()
                        if len(content) > 9:
                            sample_names = content[9:]
                    continue

                if not len(line.strip()):
                    continue

                variant = readVCF.readVCFLine(line)
                if args.passonly:
                    FILTER = line.split("\t")[6]
                    if not (FILTER in ["PASS", "."]):
                        continue

                chrA = variant.chrA
                posA = variant.posA
                chrB = variant.chrB
                posB = variant.posB
                event_type = variant.event_type
                INFO = variant.info
                FORMAT = variant.fmt

                ci_A_lower = 0
                ci_A_upper = 0
                ci_B_lower = 0
                ci_B_upper = 0
                if "CIPOS" in INFO:
                    ci = INFO["CIPOS"].replace('(','').replace(')','').split(",")
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
                    ci = INFO["CIEND"].replace('(','').replace(')','').split(",")
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
            db.insert_many(var)

    db.create_index(name='SV', columns='(var, chrA, chrB, posA, posA, posB, posB)')
    db.create_index(name='IDX', columns='(idx)')
    db.create_index(name='CHR', columns='(chrA, chrB)')
    return sample_IDs


def main(args):
    args.db = args.prefix
    if not args.files and args.folder:
        args.files = glob.glob(os.path.join(args.folder, "*.vcf")) + glob.glob(os.path.join(args.folder, "*.vcf.gz"))
    populate_db(args)
