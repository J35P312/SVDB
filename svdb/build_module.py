import glob
import logging
import os
import zlib
from pathlib import Path

from . import database, read_vcf, vcf_utils

logger = logging.getLogger(__name__)


def _store_ins_seq(seq, max_len):
    """Cap (if max_len set) then compress ins_seq for storage as a BLOB.

    Returns None when seq is absent or capped, so SVLEN-ratio matching still
    applies.  Callers that read the value back must call
    ins_similarity.decompress_ins_seq() to recover the string.
    """
    if not seq:
        return None
    if max_len is not None and len(seq) > max_len:
        return None
    return zlib.compress(seq.encode(), 6)


def populate_db(args):
    with database.DB(args.db) as db:
        tables = db.tables

        idx = 0
        if "SVDB" not in tables:
            db.create(database.CREATE_TABLE_SQL)
            sample_IDs = []
        else:
            db.drop("DROP INDEX IF EXISTS SV")
            db.drop("DROP INDEX IF EXISTS SV_POS")
            db.drop("DROP INDEX IF EXISTS IDX")
            db.drop("DROP INDEX IF EXISTS CHR")

            sample_IDs = db.sample_ids
            if sample_IDs:
                idx = 1 + int(db.query("SELECT MAX(idx) FROM SVDB")[0][0])

        db.create_ins_table()
        max_ins_seq_len = getattr(args, 'max_ins_seq_len', None)

        # populate the tables
        for vcf in args.files:
            sample_name = Path(vcf).stem.replace(".", "_")
            sample_IDs.append(sample_name)
            A = f'SELECT sample FROM SVDB WHERE sample == \'{sample_name}\' '
            hits = [hit for hit in db.query(A)]
            if hits:
                logger.debug("sample %s already in database — skipping", sample_name)
                continue
            if not os.path.exists(vcf):
                logger.warning("unable to open %s — skipping", vcf)
                continue

            var = []
            ins = []
            sample_names = []

            with vcf_utils.open_vcf(vcf) as lines:
                for line in lines:
                    if line.startswith("#"):
                        if "CHROM" in line:
                            content = line.strip().split()
                            if len(content) > 9:
                                sample_names = content[9:]
                        continue

                    if not len(line.strip()):
                        continue

                    variant = read_vcf.readVCFLine(line)
                    if args.pass_only and variant.vcf_filter not in ("PASS", "."):
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
                        ci_A_lower, ci_A_upper = vcf_utils.parse_ci(INFO["CIPOS"])
                        ci_B_lower, ci_B_upper = ci_A_lower, ci_A_upper

                    if "CIEND" in INFO:
                        ci_B_lower, ci_B_upper = vcf_utils.parse_ci(INFO["CIEND"])

                    is_ins = "INS" in event_type

                    if "GT" not in FORMAT or not len(sample_names):
                        var.append((event_type, chrA, chrB, posA, ci_A_lower,
                                    ci_A_upper, posB, ci_B_lower, ci_B_upper, sample_name, idx))
                        if is_ins:
                            ins.append((idx, _store_ins_seq(variant.ins_seq, max_ins_seq_len), variant.svlen))
                        idx += 1
                    else:
                        sample_index = 0
                        for genotype in FORMAT["GT"]:
                            if genotype not in ["0/0", "./."]:
                                var.append((event_type, chrA, chrB, posA, ci_A_lower, ci_A_upper,
                                            posB, ci_B_lower, ci_B_upper, sample_names[sample_index], idx))
                                if is_ins:
                                    ins.append((idx, _store_ins_seq(variant.ins_seq, max_ins_seq_len), variant.svlen))
                                idx += 1
                            sample_index += 1

            if var:
                db.insert_many(var)
            if ins:
                db.insert_ins_many(ins)

        db.create_index(name='SV', columns='(var, chrA, chrB, posA)')
        db.create_index(name='SV_POS', columns='(posA, posB)')
        db.create_index(name='IDX', columns='(idx)')
        db.create_index(name='CHR', columns='(chrA, chrB)')
    return sample_IDs


def upgrade_db(args):
    """Create the INS table and backfill from provided VCFs if given."""
    with database.DB(args.db) as db:
        if "SVDB" not in db.tables:
            logger.error("no SVDB table found in %s — build a database first", args.db)
            return

        if db.has_ins_table():
            logger.info("database schema is up to date")
            return

        db.create_ins_table()
        logger.info("INS table created")

        if not args.files and not args.folder:
            return

        # Backfill INS table from the provided VCFs by matching idx via SVDB lookup
        max_ins_seq_len = getattr(args, 'max_ins_seq_len', None)
        ins = []
        for vcf in args.files:
            sample_name = Path(vcf).stem.replace(".", "_")
            if not os.path.exists(vcf):
                logger.warning("unable to open %s — skipping", vcf)
                continue

            with vcf_utils.open_vcf(vcf) as lines:
                for line in lines:
                    if line.startswith("#"):
                        continue
                    if not line.strip():
                        continue

                    variant = read_vcf.readVCFLine(line)
                    if "INS" not in variant.event_type:
                        continue

                    rows = db.query(
                        f"SELECT idx FROM SVDB WHERE sample == '{sample_name}' "
                        f"AND var == '{variant.event_type}' "
                        f"AND chrA == '{variant.chrA}' AND posA == {variant.posA}"
                    )
                    for (idx,) in rows:
                        ins.append((idx, _store_ins_seq(variant.ins_seq, max_ins_seq_len), variant.svlen))

        if ins:
            db.insert_ins_many(ins)
            logger.info("backfilled %d INS entries", len(ins))


def main(args):
    args.db = args.prefix

    if getattr(args, "upgrade", False):
        if not args.files and args.folder:
            args.files = glob.glob(os.path.join(args.folder, "*.vcf")) + glob.glob(os.path.join(args.folder, "*.vcf.gz"))
        upgrade_db(args)
        return

    if not args.files and args.folder:
        args.files = glob.glob(os.path.join(args.folder, "*.vcf")) + glob.glob(os.path.join(args.folder, "*.vcf.gz"))
        if not args.files:
            logger.warning("no VCF files found in folder: %s", args.folder)
    populate_db(args)
