"""End-to-end tests for the SQLite INS table feature.

Covers: build populates INS table, --upgrade on old/current DBs, query --sqdb
warning + sequence gates, export ALT column with sequence vs symbolic fallback.
"""

import subprocess
import sys
from pathlib import Path

from svdb.database import DB, CREATE_TABLE_SQL
from svdb.ins_similarity import decompress_ins_seq

FIXTURES = Path(__file__).parent / "fixtures"
DEL_VCF = FIXTURES / "manta_chr1_del.vcf"
SVDB = [sys.executable, "-m", "svdb"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_ins_vcf(path: Path, variants: list) -> None:
    """Write an INFO-only INS VCF (no sample columns).

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


def make_gt_vcf(path: Path, sample_name: str, variants: list) -> None:
    """Write an INS VCF with a named GT sample column.

    variants: list of (chrom, pos, id, sequence) tuples.
    The sample column name differs from the filename stem to let tests verify
    that build/upgrade uses the header name, not the stem.
    """
    lines = [
        "##fileformat=VCFv4.1",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}",
    ]
    for chrom, pos, id_, seq in variants:
        lines.append(
            f"{chrom}\t{pos}\t{id_}\tN\tN{seq}\t.\tPASS\t"
            f"SVTYPE=INS;SVLEN={len(seq)};END={pos}\tGT\t0/1"
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
        assert decompress_ins_seq(rows[0][0]) == "ACGTACGT"
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
        assert decompress_ins_seq(short_seq) == "ACGT" * 10
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

    def test_upgrade_without_files_errors(self, tmp_path):
        prefix = tmp_path / "svdb"
        make_old_db(prefix, [("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, "sample_A", 0)])
        r = run("--build", "--upgrade", "--prefix", str(prefix))
        assert r.returncode != 0
        assert "requires --files or --folder" in r.stderr

    def test_upgrade_with_files_creates_ins_table_on_old_db(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGT")])
        prefix = tmp_path / "svdb"
        make_old_db(prefix, [("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, "sample_A", 0)])
        with DB(str(prefix)) as db:
            assert not db.has_ins_table()

        r = run("--build", "--upgrade", "--files", str(vcf), "--prefix", str(prefix))
        assert r.returncode == 0
        assert "INS table created" in r.stderr
        with DB(str(prefix)) as db:
            assert db.has_ins_table()

    def test_upgrade_idempotent_when_already_current(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGT")])
        build(tmp_path / "svdb", vcf)

        r = run("--build", "--upgrade", "--files", str(vcf), "--prefix", str(tmp_path / "svdb"))
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
        assert decompress_ins_seq(rows[0][0]) == "ACGTACGT"
        assert rows[0][1] == 8

    def test_upgrade_on_missing_db_exits_with_error(self, tmp_path):
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [("1", 1000, "ins1", "ACGT")])
        r = run("--build", "--upgrade", "--files", str(vcf), "--prefix", str(tmp_path / "nonexistent"))
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

    def test_query_position_only_skips_sequence_gate(self, tmp_path):
        """--data_profile position_only should match on position only, ignoring dissimilar sequences."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "AAAAAAAAAAAAAAAA")])
        make_ins_vcf(query_vcf, [("1", 1000, "q_ins", "CCCCCCCCCCCCCCCC")])
        build(tmp_path / "svdb", db_vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                "--query_vcf", str(query_vcf), "--data_profile", "position_only")
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

    def test_query_max_ins_seq_len_bypasses_sequence_gate(self, tmp_path):
        """--max_ins_seq_len caps long sequences to empty, skipping the gate → hit counted."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        # completely different 200 bp sequences → default sequence gate blocks the match
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "A" * 200)])
        make_ins_vcf(query_vcf, [("1", 1000, "q_ins", "C" * 200)])
        build(tmp_path / "svdb", db_vcf)
        # without cap: sim≈0 < 0.75 → no match
        r_nocap = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                      "--query_vcf", str(query_vcf))
        assert r_nocap.returncode == 0
        assert "OCC=" not in data_lines(r_nocap.stdout)[0]
        # with cap below 200 bp: both sequences capped to "" → gate skipped → match
        r_cap = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                    "--query_vcf", str(query_vcf), "--max_ins_seq_len", "100")
        assert r_cap.returncode == 0
        assert "OCC=1" in data_lines(r_cap.stdout)[0]

    def test_query_max_ins_seq_len_above_length_still_compares(self, tmp_path):
        """Cap above the actual sequence length → sequences still compared → gate blocks."""
        db_vcf = tmp_path / "db.vcf"
        query_vcf = tmp_path / "query.vcf"
        make_ins_vcf(db_vcf, [("1", 1000, "db_ins", "A" * 50)])
        make_ins_vcf(query_vcf, [("1", 1000, "q_ins", "C" * 50)])
        build(tmp_path / "svdb", db_vcf)
        r = run("--query", "--sqdb", str(tmp_path / "svdb.db"),
                "--query_vcf", str(query_vcf), "--max_ins_seq_len", "200")
        assert r.returncode == 0
        assert "OCC=" not in data_lines(r.stdout)[0]

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
                "--prefix", str(tmp_path / "out"), "--data_profile", "position_only")
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        # All four cluster together (position_only); most common seq is ACGTACGT
        assert any("NACGTACGT" in line for line in lines)

    def test_export_n_heavy_sequence_treated_as_symbolic(self, tmp_path):
        """N-heavy sequence is nulled out for comparison (issue #95): a clean
        sequence and an N-heavy one at the same position/SVLEN merge on
        position+SVLEN alone -- without the guard, the sequence gate would
        reject them (a clean sequence vs all-N looks maximally dissimilar)
        and they'd export as two separate lines instead of one mixed,
        symbolic cluster.
        """
        vcf = tmp_path / "sample.vcf"
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", "A" * 16),
            ("1", 1001, "ins2", "N" * 16),
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"))
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        ins_lines = [ln for ln in lines if "SVTYPE=INS" in ln]
        assert len(ins_lines) == 1
        assert "<INS>" in ins_lines[0]

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
                "--prefix", str(tmp_path / "out"), "--data_profile", "position_only")
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
                "--data_profile", "position_only", "--ins_svlen_ratio", "0.10")
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

    def test_export_max_ins_seq_len_bypasses_sequence_gate(self, tmp_path):
        """Sequences above --max_ins_seq_len are not compared: dissimilar long sequences
        at nearby positions merge on position+SVLEN only."""
        vcf = tmp_path / "sample.vcf"
        seq_a = "ACGT" * 50   # 200 bp
        seq_b = "TTTT" * 50   # 200 bp — completely different (sim ≈ 0)
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", seq_a),
            ("1", 1001, "ins2", seq_b),
        ])
        build(tmp_path / "svdb", vcf)
        # without cap: sequence gate keeps them separate
        run("--export", "--db", str(tmp_path / "svdb.db"),
            "--prefix", str(tmp_path / "out_nocap"))
        assert len(data_lines((tmp_path / "out_nocap.vcf").read_text())) == 2
        # with cap below sequence length: gate bypassed → merge into one cluster
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out_cap"),
                "--max_ins_seq_len", "100")
        assert r.returncode == 0
        assert len(data_lines((tmp_path / "out_cap.vcf").read_text())) == 1

    def test_export_max_ins_seq_len_alt_falls_back_to_symbolic(self, tmp_path):
        """Capped sequences produce <INS> in ALT (not full sequence); SVLEN preserved."""
        vcf = tmp_path / "sample.vcf"
        seq = "ACGT" * 50   # 200 bp
        make_ins_vcf(vcf, [
            ("1", 1000, "ins1", seq),
            ("1", 1001, "ins2", seq),
        ])
        build(tmp_path / "svdb", vcf)
        r = run("--export", "--db", str(tmp_path / "svdb.db"),
                "--prefix", str(tmp_path / "out"),
                "--max_ins_seq_len", "100")
        assert r.returncode == 0
        lines = data_lines((tmp_path / "out.vcf").read_text())
        assert len(lines) == 1
        fields = lines[0].split("\t")
        assert fields[4] == "<INS>"
        assert "SVLEN=200" in fields[7]

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


# ---------------------------------------------------------------------------
# Sample name consistency: build and upgrade use header names, not filename stems
# ---------------------------------------------------------------------------

# Five distinct sample names and sequences for fixture construction.
_GT_SAMPLES = [
    ("patient_A", "ACGTACGT"),
    ("patient_B", "TGCATGCA"),
    ("patient_C", "AAAACCCC"),
    ("patient_D", "GGGGTTTT"),
    ("patient_E", "ATATATAT"),
]


class TestBuildSampleNames:

    def test_build_uses_header_name_not_stem(self, tmp_path):
        """Build stores the header column name, not the filename stem."""
        vcf = tmp_path / "file_with_different_stem.vcf"
        make_gt_vcf(vcf, "MySample", [("1", 1000, "ins1", "ACGTACGT")])
        build(tmp_path / "svdb", vcf)
        with DB(str(tmp_path / "svdb")) as db:
            stored = db.sample_ids
        assert "MySample" in stored
        assert "file_with_different_stem" not in stored

    def test_build_dedup_uses_header_name(self, tmp_path):
        """Building the same GT-column VCF twice stores data only once."""
        vcf = tmp_path / "file.vcf"
        make_gt_vcf(vcf, "PatientX", [("1", 1000, "ins1", "ACGTACGT")])
        build(tmp_path / "svdb", vcf)
        # second build on the same DB must not double-insert
        build(tmp_path / "svdb", vcf)
        with DB(str(tmp_path / "svdb")) as db:
            count = db.query("SELECT COUNT(*) FROM SVDB")[0][0]
        assert count == 1


class TestUpgradeSampleNames:

    def _build_old_db(self, prefix: Path, samples_seqs: list) -> None:
        """Create an old DB (no INS table) from (sample_name, seq) pairs at distinct positions."""
        rows = [
            ("INS", "1", "1", 1000 + i * 100, 0, 0, 1000 + i * 100, 0, 0, sname, i)
            for i, (sname, _) in enumerate(samples_seqs)
        ]
        make_old_db(prefix, rows)

    def _make_gt_vcfs(self, tmp_path: Path, samples_seqs: list) -> list:
        """Write one GT-column VCF per (sample_name, seq) pair; filename stems differ from names."""
        paths = []
        for i, (sname, seq) in enumerate(samples_seqs):
            vcf = tmp_path / f"vcf_{i:02d}.vcf"
            make_gt_vcf(vcf, sname, [("1", 1000 + i * 100, "ins1", seq)])
            paths.append(vcf)
        return paths

    def test_upgrade_uses_header_name_not_stem(self, tmp_path):
        """upgrade_db matches SVDB rows by header sample name, not filename stem."""
        sname, seq = "MySample", "ACGTACGT"
        vcf = tmp_path / "completely_different_stem.vcf"
        make_gt_vcf(vcf, sname, [("1", 1000, "ins1", seq)])
        make_old_db(tmp_path / "svdb", [
            ("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, sname, 0)
        ])
        r = run("--build", "--upgrade", "--files", str(vcf), "--prefix", str(tmp_path / "svdb"))
        assert r.returncode == 0
        with DB(str(tmp_path / "svdb")) as db:
            rows = db.query("SELECT ins_seq FROM INS")
        assert len(rows) == 1, "header-name match failed — stem was used instead"

    def test_upgrade_5_samples_all_backfilled(self, tmp_path):
        """Upgrading with all 5 VCFs backfills every sample."""
        self._build_old_db(tmp_path / "svdb", _GT_SAMPLES)
        vcf_paths = self._make_gt_vcfs(tmp_path, _GT_SAMPLES)
        r = run(
            "--build", "--upgrade",
            "--files", *[str(v) for v in vcf_paths],
            "--prefix", str(tmp_path / "svdb"),
        )
        assert r.returncode == 0
        with DB(str(tmp_path / "svdb")) as db:
            count = db.query("SELECT COUNT(*) FROM INS")[0][0]
        assert count == 5

    def test_upgrade_partial_warns_missing_samples(self, tmp_path):
        """Upgrading with only 3 of 5 VCFs warns about the 2 missing samples."""
        self._build_old_db(tmp_path / "svdb", _GT_SAMPLES)
        vcf_paths = self._make_gt_vcfs(tmp_path, _GT_SAMPLES)
        r = run(
            "--build", "--upgrade",
            "--files", *[str(v) for v in vcf_paths[:3]],
            "--prefix", str(tmp_path / "svdb"),
        )
        assert r.returncode == 0
        # 3 samples backfilled
        with DB(str(tmp_path / "svdb")) as db:
            count = db.query("SELECT COUNT(*) FROM INS")[0][0]
        assert count == 3
        # 2 missing-sample warnings
        warnings = [ln for ln in r.stderr.splitlines() if "INS data not backfilled" in ln]
        assert len(warnings) == 2

    def test_upgrade_extra_vcf_logs_no_db_entries(self, tmp_path):
        """A VCF whose sample has no SVDB rows produces an info log, not an error."""
        sname, seq = _GT_SAMPLES[0]
        vcf_known = tmp_path / "vcf_00.vcf"
        make_gt_vcf(vcf_known, sname, [("1", 1000, "ins1", seq)])
        make_old_db(tmp_path / "svdb", [
            ("INS", "1", "1", 1000, 0, 0, 1000, 0, 0, sname, 0)
        ])

        extra_vcf = tmp_path / "extra.vcf"
        make_gt_vcf(extra_vcf, "UnknownSample", [("1", 2000, "ins1", "TTTTTTTT")])

        r = run(
            "--build", "--upgrade",
            "--files", str(vcf_known), str(extra_vcf),
            "--prefix", str(tmp_path / "svdb"),
        )
        assert r.returncode == 0
        assert "UnknownSample" in r.stderr
        assert "no entries in the database" in r.stderr
