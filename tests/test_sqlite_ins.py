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
    with DB(str(prefix)) as db:
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
        with DB(str(tmp_path / "svdb")) as db:
            assert db.has_ins_table()

    def test_build_populates_ins_seq_and_len(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGTACGT")])
        build(tmp_path / "svdb", vcf)
        with DB(str(tmp_path / "svdb")) as db:
            rows = db.query("SELECT ins_seq, ins_len FROM INS")
        assert len(rows) == 1
        assert rows[0][0] == "ACGTACGT"
        assert rows[0][1] == 8

    def test_build_del_variants_not_in_ins_table(self, tmp_path):
        build(tmp_path / "svdb", DEL_VCF)
        with DB(str(tmp_path / "svdb")) as db:
            assert db.has_ins_table()
            assert db.query("SELECT COUNT(*) FROM INS")[0][0] == 0

    def test_build_max_ins_seq_len_caps_long_sequences(self, tmp_path):
        """Sequences above --max_ins_seq_len are stored as NULL; ins_len is preserved."""
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [
            ("1", 1000, "short", "ACGT" * 10),   # 40 bp — below cap
            ("1", 2000, "long",  "ACGT" * 200),  # 800 bp — above cap
        ])
        subprocess.run(
            SVDB + ["--build", "--files", str(vcf), "--prefix", str(tmp_path / "svdb"),
                    "--max_ins_seq_len", "100"],
            capture_output=True, check=True,
        )
        with DB(str(tmp_path / "svdb")) as db:
            rows = db.query("SELECT ins_seq, ins_len FROM INS ORDER BY ins_len")
        short_seq, short_len = rows[0]
        long_seq, long_len = rows[1]
        assert short_seq == "ACGT" * 10
        assert short_len == 40
        assert long_seq is None       # capped → no sequence stored
        assert long_len == 800        # SVLEN still preserved for ratio matching

    def test_build_multiple_ins_variants_all_stored(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "ACGT"),
            ("1", 2000, "ins2", "TTTT"),
            ("1", 3000, "ins3", "CCCC"),
        ])
        build(tmp_path / "svdb", vcf)
        with DB(str(tmp_path / "svdb")) as db:
            assert db.query("SELECT COUNT(*) FROM INS")[0][0] == 3


# ---------------------------------------------------------------------------
# --upgrade
# ---------------------------------------------------------------------------

class TestUpgrade:

    def test_upgrade_creates_ins_table_on_old_db(self, tmp_path):
        prefix = tmp_path / "svdb"
        make_old_db(prefix, [("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, "sample_A", 0)])
        with DB(str(prefix)) as db:
            assert not db.has_ins_table()

        r = run("--build", "--upgrade", "--prefix", str(prefix))
        assert r.returncode == 0
        assert "INS table created" in r.stderr
        with DB(str(prefix)) as db:
            assert db.has_ins_table()

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

        with DB(str(tmp_path / "svdb")) as db:
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

    def test_query_svlen_ratio_blocks_dissimilar_length(self, tmp_path):
        """SVLEN ratio gate: DB SVLEN=80, query SVLEN=100 → ratio=0.80 < default 0.90."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "A" * 80)])
        make_ins_vcf(query_vcf, [("1", 1000, "q_ins", "A" * 100)])
        build(tmp_path / "svdb", db_vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                "--query_vcf", str(query_vcf))
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert len(lines) == 1
        assert "OCC=" not in lines[0]

    def test_query_svlen_ratio_permissive_allows_match(self, tmp_path):
        """Lowering --ins_svlen_ratio lets the same pair through."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "A" * 80)])
        make_ins_vcf(query_vcf, [("1", 1000, "q_ins", "A" * 100)])
        build(tmp_path / "svdb", db_vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                "--query_vcf", str(query_vcf), "--ins_svlen_ratio", "0.75")
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert len(lines) == 1
        assert "OCC=1" in lines[0]

    def test_query_ins_distance_blocks_far_insertion(self, tmp_path):
        """Insertion 30 bp away is outside the default ins_distance=25 window."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "ACGTACGT")])
        make_ins_vcf(query_vcf, [("1", 1030, "q_ins", "ACGTACGT")])
        build(tmp_path / "svdb", db_vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                "--query_vcf", str(query_vcf))
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert len(lines) == 1
        assert "OCC=" not in lines[0]

    def test_query_ins_distance_custom_allows_far_insertion(self, tmp_path):
        """--ins_distance 50 brings a 30 bp gap within range."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "ACGTACGT")])
        make_ins_vcf(query_vcf, [("1", 1030, "q_ins", "ACGTACGT")])
        build(tmp_path / "svdb", db_vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                "--query_vcf", str(query_vcf), "--ins_distance", "50")
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert len(lines) == 1
        assert "OCC=1" in lines[0]

    def test_query_data_profile_sample_blocks_moderate_similarity(self, tmp_path):
        """--data_profile sample sets threshold=0.85; sequence with sim≈0.80 should not match."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        # 50 A's in DB; 40 A's + 10 C's in query → Levenshtein sim = 1 - 10/50 = 0.80
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "A" * 50)])
        make_ins_vcf(query_vcf, [("1", 1000, "q_ins", "A" * 40 + "C" * 10)])
        build(tmp_path / "svdb", db_vcf)
        # Default threshold 0.75: sim 0.80 passes → should match
        r_default = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                        "--query_vcf", str(query_vcf))
        assert "OCC=1" in data_lines(r_default.stdout)[0]
        # sample profile threshold 0.85: sim 0.80 fails → should not match
        r_profile = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                        "--query_vcf", str(query_vcf), "--data_profile", "sample")
        assert r_profile.returncode == 0
        assert "OCC=" not in data_lines(r_profile.stdout)[0]


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

    def test_export_ins_svlen_ratio_separates_different_lengths(self, tmp_path):
        """Default ratio=0.90 keeps very different lengths in separate clusters."""
        vcf = tmp_path / "sample.vcf"
        # ratio = 20/100 = 0.20, well below default 0.90
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "A" * 20),
            ("1", 1001, "ins2", "A" * 100),
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"), "--no_ins_seq")
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        assert len(lines) == 2

    def test_export_ins_svlen_ratio_permissive_merges_different_lengths(self, tmp_path):
        """Lowering --ins_svlen_ratio lets differently-sized insertions cluster."""
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "A" * 20),
            ("1", 1001, "ins2", "A" * 100),
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"),
                "--no_ins_seq", "--ins_svlen_ratio", "0.10")
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        assert len(lines) == 1

    def test_export_ins_distance_separates_far_insertions(self, tmp_path):
        """Two insertions 30 bp apart stay separate with the default ins_distance=25."""
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "ACGTACGT"),
            ("1", 1030, "ins2", "ACGTACGT"),
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"))
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        assert len(lines) == 2

    def test_export_ins_distance_custom_merges_far_insertions(self, tmp_path):
        """--ins_distance 50 brings two insertions 30 bp apart into the same cluster."""
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "ACGTACGT"),
            ("1", 1030, "ins2", "ACGTACGT"),
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"), "--ins_distance", "50")
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        assert len(lines) == 1

    def test_export_data_profile_sample_separates_moderate_similarity(self, tmp_path):
        """--data_profile sample (threshold=0.85) keeps sim≈0.80 insertions separate."""
        vcf = tmp_path / "sample.vcf"
        # 50 A's vs 40 A's + 10 C's → Levenshtein sim = 1 - 10/50 = 0.80
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "A" * 50),
            ("1", 1001, "ins2", "A" * 40 + "C" * 10),
        ])
        build(tmp_path / "svdb", vcf)
        # Default threshold 0.75: sim 0.80 passes → one cluster
        run("--export", "--db", str(tmp_path / "svdb.db"),
            "--prefix", str(tmp_path / "out_default"))
        assert len(data_lines((tmp_path / "out_default.vcf").read_text())) == 1
        # sample profile threshold 0.85: sim 0.80 fails → two clusters
        r_profile = run("--export", "--db", str(tmp_path / "svdb.db"),
                        "--prefix", str(tmp_path / "out_profile"),
                        "--data_profile", "sample")
        assert r_profile.returncode == 0
        assert len(data_lines((tmp_path / "out_profile.vcf").read_text())) == 2

    def test_export_subclusters_get_distinct_ids(self, tmp_path):
        """Each sub-cluster produced by overlap_cluster must get a unique cluster ID.

        Regression for a pre-existing bug where i was not incremented inside the
        overlap_cluster loop — all sub-clusters from one DBSCAN group shared the
        same cluster_N id.
        """
        vcf = tmp_path / "sample.vcf"
        # Two insertions within ins_distance (5 bp apart) but with completely
        # different sequences (similarity ≈ 0) → sequence gate splits them.
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "ACGT" * 10),   # 40 bp, all ACGT
            ("1", 1005, "ins2", "TTTT" * 10),   # 40 bp, all T — sim ≈ 0.25
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"))
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        assert len(lines) == 2
        ids = [line.split("\t")[2] for line in lines]
        assert ids[0] != ids[1], f"both sub-clusters got the same ID: {ids[0]}"
