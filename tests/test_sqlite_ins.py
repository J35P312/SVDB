"""End-to-end tests for the SQLite INS table feature.

Covers: build populates INS table, --upgrade on old/current DBs, query --sqdb
warning + sequence gates, export ALT column with sequence vs symbolic fallback.
"""

import subprocess
import sys
from pathlib import Path

from svdb.database import DB, CREATE_TABLE_SQL

FIXTURES = Path(__file__).parent / "fixtures"
DEL_VCF = FIXTURES / "manta_chr1_del.vcf"
SVDB = [sys.executable, "-m", "svdb"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_ins_vcf(path: Path, variants: list) -> None:
    """Write a minimal single-sample INS VCF.

    variants: list of (chrom, pos, id, sequence) tuples.
    """
    lines = [
        "##fileformat=VCFv4.1",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for chrom, pos, id_, seq in variants:
        lines.append(
            f"{chrom}\t{pos}\t{id_}\tN\tN{seq}\t.\tPASS\t"
            f"SVTYPE=INS;SVLEN={len(seq)};END={pos}"
        )
    path.write_text("\n".join(lines) + "\n")


def make_old_db(prefix: Path, svdb_rows: list) -> None:
    """Create a database with only the SVDB table (no INS table) — simulates pre-upgrade DB."""
    db = DB(str(prefix))
    db.create(CREATE_TABLE_SQL)
    db.insert_many(svdb_rows)


def build(prefix: Path, *vcfs: Path) -> None:
    subprocess.run(
        SVDB + ["--build", "--files"] + [str(v) for v in vcfs] + ["--prefix", str(prefix)],
        capture_output=True, check=True,
    )


def run(*args) -> subprocess.CompletedProcess:
    return subprocess.run(SVDB + list(args), capture_output=True, text=True)


def data_lines(text: str) -> list:
    return [line for line in text.splitlines() if line and not line.startswith("#")]


# ---------------------------------------------------------------------------
# Build: INS table population
# ---------------------------------------------------------------------------

class TestBuildINSTable:

    def test_build_creates_ins_table(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGTACGT")])
        build(tmp_path / "svdb", vcf)
        assert DB(str(tmp_path / "svdb")).has_ins_table()

    def test_build_populates_ins_seq_and_len(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGTACGT")])
        build(tmp_path / "svdb", vcf)
        db = DB(str(tmp_path / "svdb"))
        rows = db.query("SELECT ins_seq, ins_len FROM INS")
        assert len(rows) == 1
        assert rows[0][0] == "ACGTACGT"
        assert rows[0][1] == 8

    def test_build_del_variants_not_in_ins_table(self, tmp_path):
        build(tmp_path / "svdb", DEL_VCF)
        db = DB(str(tmp_path / "svdb"))
        assert db.has_ins_table()
        assert db.query("SELECT COUNT(*) FROM INS")[0][0] == 0

    def test_build_multiple_ins_variants_all_stored(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "ACGT"),
            ("1", 2000, "ins2", "TTTT"),
            ("1", 3000, "ins3", "CCCC"),
        ])
        build(tmp_path / "svdb", vcf)
        db = DB(str(tmp_path / "svdb"))
        assert db.query("SELECT COUNT(*) FROM INS")[0][0] == 3


# ---------------------------------------------------------------------------
# --upgrade
# ---------------------------------------------------------------------------

class TestUpgrade:

    def test_upgrade_creates_ins_table_on_old_db(self, tmp_path):
        prefix = tmp_path / "svdb"
        make_old_db(prefix, [("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, "sample_A", 0)])
        assert not DB(str(prefix)).has_ins_table()

        r = run("--build", "--upgrade", "--prefix", str(prefix))
        assert r.returncode == 0
        assert "INS table created" in r.stderr
        assert DB(str(prefix)).has_ins_table()

    def test_upgrade_idempotent_when_already_current(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGT")])
        build(tmp_path / "svdb", vcf)

        r = run("--build", "--upgrade", "--prefix", str(tmp_path / "svdb"))
        assert r.returncode == 0
        assert "up to date" in r.stderr

    def test_upgrade_with_files_backfills_ins_data(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGTACGT")])
        # build old DB without INS table, with one matching SVDB row
        make_old_db(tmp_path / "svdb", [
            ("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, "sample", 0)
        ])

        r = run("--build", "--upgrade", "--files", str(vcf), "--prefix", str(tmp_path / "svdb"))
        assert r.returncode == 0

        db = DB(str(tmp_path / "svdb"))
        rows = db.query("SELECT ins_seq, ins_len FROM INS")
        assert len(rows) == 1
        assert rows[0][0] == "ACGTACGT"
        assert rows[0][1] == 8

    def test_upgrade_on_missing_db_exits_gracefully(self, tmp_path):
        r = run("--build", "--upgrade", "--prefix", str(tmp_path / "nonexistent"))
        assert r.returncode == 0
        assert "no SVDB table" in r.stderr


# ---------------------------------------------------------------------------
# Query --sqdb: warning + sequence gates
# ---------------------------------------------------------------------------

class TestQuerySQLiteINS:

    def test_query_warns_when_ins_table_absent(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGT")])
        make_old_db(tmp_path / "svdb", [
            ("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, "sample_A", 0)
        ])
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"), "--query_vcf", str(vcf))
        assert r.returncode == 0
        assert "does not contain insertion sequence" in r.stderr

    def test_query_no_warning_when_ins_table_present(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGT")])
        build(tmp_path / "svdb", vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"), "--query_vcf", str(vcf))
        assert r.returncode == 0
        assert "does not contain insertion sequence" not in r.stderr

    def test_query_counts_matching_insertion(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGTACGT")])
        build(tmp_path / "svdb", vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"), "--query_vcf", str(vcf))
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert len(lines) == 1
        assert "OCC=1" in lines[0]

    def test_query_sequence_gate_blocks_dissimilar_insertion(self, tmp_path):
        """Completely different sequence at same position should not be counted."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "AAAAAAAAAAAAAAAA")])
        make_ins_vcf(query_vcf, [("1", 1000, "q_ins", "CCCCCCCCCCCCCCCC")])
        build(tmp_path / "svdb", db_vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                "--query_vcf", str(query_vcf), "--ins_seq_similarity", "0.75")
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert len(lines) == 1
        assert "OCC=" not in lines[0]

    def test_query_no_ins_seq_skips_sequence_gate(self, tmp_path):
        """--no_ins_seq should match on position only, ignoring dissimilar sequences."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "AAAAAAAAAAAAAAAA")])
        make_ins_vcf(query_vcf, [("1", 1000, "q_ins", "CCCCCCCCCCCCCCCC")])
        build(tmp_path / "svdb", db_vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                "--query_vcf", str(query_vcf), "--no_ins_seq")
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert len(lines) == 1
        assert "OCC=1" in lines[0]


# ---------------------------------------------------------------------------
# Export: ALT column + warning
# ---------------------------------------------------------------------------

class TestExportINS:

    def test_export_alt_uses_sequence_when_ins_table_present(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        # two identical insertions to form a cluster
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "ACGTACGT"),
            ("1", 1001, "ins2", "ACGTACGT"),
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"))
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        assert any("NACGTACGT" in line for line in lines)

    def test_export_alt_falls_back_to_symbolic_without_ins_table(self, tmp_path):
        make_old_db(tmp_path / "svdb", [
            ("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, "sample_A", 0)
        ])
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"))
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        assert any("<INS>" in line for line in lines)

    def test_export_warns_when_ins_table_absent(self, tmp_path):
        make_old_db(tmp_path / "svdb", [
            ("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, "sample_A", 0)
        ])
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"))
        assert r.returncode == 0
        assert "does not contain insertion sequence" in r.stderr

    def test_export_no_warning_when_ins_table_present(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGT")])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"))
        assert r.returncode == 0
        assert "does not contain insertion sequence" not in r.stderr

    def test_export_cluster_representative_uses_most_common_sequence(self, tmp_path):
        """Three identical sequences vs one different — most common should win."""
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "ACGTACGT"),
            ("1", 1001, "ins2", "ACGTACGT"),
            ("1", 1002, "ins3", "ACGTACGT"),
            ("1", 1003, "ins4", "TTTTTTTT"),
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"), "--no_ins_seq")
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        # All four cluster together (--no_ins_seq); most common seq is ACGTACGT
        assert any("NACGTACGT" in line for line in lines)
