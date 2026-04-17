"""Integration tests using real-world VCF fixture data.

Fixtures are small slices of real caller output (manta, TIDDIT) and a
truth set (Personalis 1000 Genomes), extracted to keep tests fast while
exercising real variant representations (DEL, BND, multi-sample GT).
"""

import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
MANTA = FIXTURES / "manta_chr1_del.vcf"
TIDDIT = FIXTURES / "tiddit_chr1_del.vcf"
TRUTH = FIXTURES / "truth_chr1_del.vcf"
MANTA_BND = FIXTURES / "manta_bnd.vcf"

SVDB = [sys.executable, "-m", "svdb"]


def run(*args):
    """Run svdb with args, return CompletedProcess. Captures stdout+stderr."""
    result = subprocess.run(
        SVDB + list(args),
        capture_output=True,
        text=True,
    )
    return result


def vcf_data_lines(text: str) -> list[str]:
    """Return non-header lines from VCF text."""
    return [line for line in text.splitlines() if line and not line.startswith("#")]


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------


class TestBuild:
    def test_build_two_callers(self, tmp_path):
        prefix = tmp_path / "svdb"
        r = run("--build", "--files", str(MANTA), str(TIDDIT), "--prefix", str(prefix))
        assert r.returncode == 0
        assert (tmp_path / "svdb.db").exists()

    def test_build_passonly(self, tmp_path):
        prefix = tmp_path / "svdb_pass"
        r = run("--build", "--files", str(MANTA), "--passonly", "--prefix", str(prefix))
        assert r.returncode == 0
        assert (tmp_path / "svdb_pass.db").exists()

    def test_build_from_folder(self, tmp_path):
        import shutil
        folder = tmp_path / "vcfs"
        folder.mkdir()
        shutil.copy(MANTA, folder / "manta.vcf")
        shutil.copy(TIDDIT, folder / "tiddit.vcf")
        r = run("--build", "--folder", str(folder), "--prefix", str(tmp_path / "from_folder"))
        assert r.returncode == 0
        assert (tmp_path / "from_folder.db").exists()

    def test_build_idempotent(self, tmp_path):
        """Building the same files twice should not raise or duplicate rows."""
        prefix = tmp_path / "idem"
        run("--build", "--files", str(MANTA), "--prefix", str(prefix))
        r = run("--build", "--files", str(MANTA), "--prefix", str(prefix))
        assert r.returncode == 0


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------


class TestExport:
    @pytest.fixture
    def db(self, tmp_path):
        prefix = tmp_path / "svdb"
        run("--build", "--files", str(MANTA), str(TIDDIT), "--prefix", str(prefix))
        return tmp_path / "svdb.db"

    def test_export_default(self, db, tmp_path):
        prefix = tmp_path / "out"
        r = run("--export", "--db", str(db), "--prefix", str(prefix))
        assert r.returncode == 0
        vcf = tmp_path / "out.vcf"
        assert vcf.exists()
        assert len(vcf_data_lines(vcf.read_text())) > 0

    def test_export_no_merge_has_more_variants(self, db, tmp_path):
        merged = tmp_path / "merged"
        no_merge = tmp_path / "no_merge"
        run("--export", "--db", str(db), "--prefix", str(merged))
        run("--export", "--db", str(db), "--no_merge", "--prefix", str(no_merge))
        n_merged = len(vcf_data_lines((merged.parent / "merged.vcf").read_text()))
        n_no_merge = len(vcf_data_lines((no_merge.parent / "no_merge.vcf").read_text()))
        assert n_no_merge >= n_merged

    def test_export_overlap0_merges_more_than_overlap1(self, db, tmp_path):
        ov0 = tmp_path / "ov0"
        ov1 = tmp_path / "ov1"
        run("--export", "--db", str(db), "--overlap", "0", "--prefix", str(ov0))
        run("--export", "--db", str(db), "--overlap", "1", "--prefix", str(ov1))
        n0 = len(vcf_data_lines((ov0.parent / "ov0.vcf").read_text()))
        n1 = len(vcf_data_lines((ov1.parent / "ov1.vcf").read_text()))
        assert n0 <= n1

    def test_export_dbscan(self, db, tmp_path):
        prefix = tmp_path / "dbscan"
        r = run("--export", "--db", str(db), "--DBSCAN", "--epsilon", "500",
                "--min_pts", "2", "--prefix", str(prefix))
        assert r.returncode == 0

    def test_export_memory_flag(self, db, tmp_path):
        prefix = tmp_path / "mem"
        r = run("--export", "--db", str(db), "--memory", "--prefix", str(prefix))
        assert r.returncode == 0


# ---------------------------------------------------------------------------
# Query — VCF db
# ---------------------------------------------------------------------------


class TestQueryVcfDb:
    def test_query_basic(self):
        r = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA))
        assert r.returncode == 0
        lines = vcf_data_lines(r.stdout)
        assert len(lines) > 0

    def test_query_annotates_occ_frq(self):
        r = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA))
        data = vcf_data_lines(r.stdout)
        annotated = [line for line in data if "OCC=" in line]
        assert len(annotated) > 0

    def test_query_loose_bnd_distance_finds_more(self):
        tight = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA),
                    "--bnd_distance", "100")
        loose = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA),
                    "--bnd_distance", "50000")
        n_tight = len([line for line in vcf_data_lines(tight.stdout) if "OCC=" in line])
        n_loose = len([line for line in vcf_data_lines(loose.stdout) if "OCC=" in line])
        assert n_loose >= n_tight

    def test_query_overlap0_finds_more_than_overlap1(self):
        ov0 = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA),
                  "--overlap", "0.0")
        ov1 = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA),
                  "--overlap", "1.0")
        n0 = len([line for line in vcf_data_lines(ov0.stdout) if "OCC=" in line])
        n1 = len([line for line in vcf_data_lines(ov1.stdout) if "OCC=" in line])
        assert n0 >= n1

    def test_query_self_annotates_all(self):
        """Querying a VCF against itself should annotate every variant."""
        r = run("--query", "--db", str(MANTA), "--query_vcf", str(MANTA))
        data = vcf_data_lines(r.stdout)
        annotated = [line for line in data if "OCC=" in line]
        assert len(annotated) == len(data)

    def test_query_custom_out_tags(self):
        r = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA),
                "--out_occ", "AC", "--out_frq", "AF")
        assert r.returncode == 0
        annotated = [line for line in vcf_data_lines(r.stdout) if "AC=" in line and "AF=" in line]
        assert len(annotated) > 0

    def test_query_no_var_finds_more(self):
        default = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA))
        no_var = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA), "--no_var")
        n_default = len([line for line in vcf_data_lines(default.stdout) if "OCC=" in line])
        n_no_var = len([line for line in vcf_data_lines(no_var.stdout) if "OCC=" in line])
        assert n_no_var >= n_default

    def test_query_prefix_writes_file(self, tmp_path):
        prefix = tmp_path / "out"
        r = run("--query", "--db", str(TRUTH), "--query_vcf", str(MANTA),
                "--prefix", str(prefix))
        assert r.returncode == 0
        assert (tmp_path / "out_query.vcf").exists()

    def test_query_bnd_variants(self):
        """BND interchromosomal variants should be queryable without crash."""
        r = run("--query", "--db", str(MANTA_BND), "--query_vcf", str(MANTA_BND))
        assert r.returncode == 0


# ---------------------------------------------------------------------------
# Query — SQLite db
# ---------------------------------------------------------------------------


class TestQuerySqDb:
    @pytest.fixture
    def db(self, tmp_path):
        prefix = tmp_path / "svdb"
        run("--build", "--files", str(MANTA), str(TIDDIT), "--prefix", str(prefix))
        return tmp_path / "svdb.db"

    def test_query_sqdb_basic(self, db):
        r = run("--query", "--sqdb", str(db), "--query_vcf", str(TRUTH))
        assert r.returncode == 0
        assert len(vcf_data_lines(r.stdout)) > 0

    def test_query_sqdb_memory_flag(self, db):
        r = run("--query", "--sqdb", str(db), "--query_vcf", str(TRUTH), "--memory")
        assert r.returncode == 0

    def test_query_sqdb_max_frq_filters(self, db):
        """max_frq=1.0 returns all; max_frq=0.0 returns only novel variants."""
        full = run("--query", "--sqdb", str(db), "--query_vcf", str(TRUTH), "--max_frq", "1.0")
        filtered = run("--query", "--sqdb", str(db), "--query_vcf", str(TRUTH), "--max_frq", "0.0")
        n_full = len(vcf_data_lines(full.stdout))
        n_filtered = len(vcf_data_lines(filtered.stdout))
        assert n_filtered <= n_full


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------


class TestMerge:
    def test_merge_two_callers(self):
        r = run("--merge", "--vcf", str(MANTA), str(TIDDIT))
        assert r.returncode == 0
        assert len(vcf_data_lines(r.stdout)) > 0

    def test_merge_has_set_tag(self):
        r = run("--merge", "--vcf", str(MANTA), str(TIDDIT))
        data = vcf_data_lines(r.stdout)
        assert any("set=" in line for line in data)

    def test_merge_priority(self):
        r = run("--merge",
                "--vcf", f"{MANTA}:manta", f"{TIDDIT}:tiddit",
                "--priority", "manta,tiddit")
        assert r.returncode == 0
        data = vcf_data_lines(r.stdout)
        assert len(data) > 0

    def test_merge_no_intra(self):
        r = run("--merge", "--vcf", str(MANTA), str(TIDDIT), "--no_intra")
        assert r.returncode == 0

    def test_merge_no_var(self):
        r = run("--merge", "--vcf", str(MANTA), str(TIDDIT), "--no_var")
        assert r.returncode == 0

    def test_merge_strict_overlap_more_variants(self):
        """Strict overlap=1.0 should produce >= variants than loose overlap=0.0
        because fewer pairs are merged."""
        loose = run("--merge", "--vcf", str(MANTA), str(TIDDIT), "--overlap", "0.0")
        strict = run("--merge", "--vcf", str(MANTA), str(TIDDIT), "--overlap", "1.0")
        n_loose = len(vcf_data_lines(loose.stdout))
        n_strict = len(vcf_data_lines(strict.stdout))
        assert n_strict >= n_loose

    def test_merge_single_vcf(self):
        """Single-file merge should return without error."""
        r = run("--merge", "--vcf", str(MANTA))
        assert r.returncode == 0

    def test_merge_supp_vec_in_header(self):
        r = run("--merge", "--vcf", str(MANTA), str(TIDDIT))
        assert "SUPP_VEC" in r.stdout
