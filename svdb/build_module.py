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


def _vcf_sample_names(vcf_path: str) -> list:
    """Return the sample names this VCF will store.

    VCFs with sample columns (len(header) > 9) use the header column names;
    VCFs without sample columns (INFO-only format) fall back to the filename stem.
    """
    stem = Path(vcf_path).stem.replace(".", "_")
    with vcf_utils.open_vcf(vcf_path) as lines:
        for line in lines:
            if not line.startswith("#"):
                break
            if "CHROM" in line:
                content = line.strip().split()
                if len(content) > 9:
                    return content[9:]
    return [stem]


def _iter_vcf_variants(vcf_path: str):
    """Yield (sample_name, variant) pairs from a VCF.

    For VCFs with sample columns and a GT FORMAT field, yields one entry per
    non-ref genotype using the corresponding header column name.  For VCFs
    without sample columns (INFO-only), yields one entry per variant line
    using the filename stem.
    """
    stem = Path(vcf_path).stem.replace(".", "_")
    sample_names: list = []

    with vcf_utils.open_vcf(vcf_path) as lines:
        for line in lines:
            if line.startswith("#"):
                if "CHROM" in line:
                    content = line.strip().split()
                    if len(content) > 9:
                        sample_names = content[9:]
                continue
            if not line.strip():
                continue

            variant = read_vcf.readVCFLine(line)
            if variant is None:
                continue

            if "GT" not in variant.fmt or not sample_names:
                yield stem, variant
            else:
                for sample_index, genotype in enumerate(variant.fmt["GT"]):
                    if genotype not in ["0/0", "./."]:
                        yield sample_names[sample_index], variant


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

        for vcf in args.files:
            if not os.path.exists(vcf):
                logger.warning("unable to open %s — skipping", vcf)
                continue

            vcf_samples = _vcf_sample_names(vcf)
            if any(db.query(f"SELECT 1 FROM SVDB WHERE sample == '{s}' LIMIT 1")
                   for s in vcf_samples):
                logger.debug("sample(s) %s already in database — skipping", vcf_samples)
                continue

            var = []
            ins = []

            for sample_name, variant in _iter_vcf_variants(vcf):
                if args.pass_only and variant.vcf_filter not in ("PASS", "."):
                    continue

                chrA = variant.chrA
                posA = variant.posA
                chrB = variant.chrB
                posB = variant.posB
                event_type = variant.event_type
                INFO = variant.info

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

                var.append((event_type, chrA, chrB, posA, ci_A_lower,
                            ci_A_upper, posB, ci_B_lower, ci_B_upper, sample_name, idx))
                if is_ins:
                    ins.append((idx, _store_ins_seq(variant.ins_seq, max_ins_seq_len), variant.svlen))
                idx += 1

            if var:
                db.insert_many(var)
                sample_IDs.extend(vcf_samples)
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

        max_ins_seq_len = getattr(args, 'max_ins_seq_len', None)
        ins = []
        db_samples = set(db.sample_ids)
        vcf_sample_names: set = set()

        for vcf in args.files:
            if not os.path.exists(vcf):
                logger.warning("unable to open %s — skipping", vcf)
                continue

            for sample_name, variant in _iter_vcf_variants(vcf):
                vcf_sample_names.add(sample_name)
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

        for s in db_samples - vcf_sample_names:
            logger.warning(
                "sample %s is in the database but no VCF was provided; INS data not backfilled", s
            )
        for s in vcf_sample_names - db_samples:
            logger.info(
                "sample %s from the provided VCF(s) has no entries in the database; skipped", s
            )


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
