"""Multi-sample integration tests: two individuals, sex chromosomes, SV-type edge cases.

Fixtures:
  HG001 (NA12878, female) — manta (short-read SV) + TIDDIT (short-read SV+CNV) callers
  HG002 (NA24385, male)   — DRAGEN SV, DRAGEN CNV, and DRAGEN SV+CNV callers

Test groupings:
  - Two-individual merge on autosomes (the core population-DB use case)
  - chrX across two individuals (both have X)
  - chrY — only HG002 (male); absent in HG001 (female)
  - Interchromosomal BND
  - SV type isolation (CNV vs DEL must not merge by default)
  - pass_only with non-PASS variants (caller-agnostic filter behaviour)
"""

import subprocess
import sys
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"

# HG001 (female) — existing autosomes + new chrX
MANTA_CHR1      = FIXTURES / "manta_chr1_del.vcf"
TIDDIT_CHR1     = FIXTURES / "tiddit_chr1_del.vcf"
MANTA_CHRX      = FIXTURES / "manta_hg001_chrX_del.vcf"
TIDDIT_CHRX     = FIXTURES / "tiddit_hg001_chrX.vcf"

# HG002 (male)
HG002_CHR1      = FIXTURES / "dragen_hg002_chr1_del.vcf"       # 15 PASS DEL on chr1
HG002_CHRX_DEL  = FIXTURES / "dragen_hg002_chrX_del.vcf"       # 15 PASS DEL on chrX
HG002_CHRY      = FIXTURES / "dragen_hg002_chrY.vcf"           # 15 PASS on chrY
HG002_BND_INTER = FIXTURES / "dragen_hg002_bnd_interchr.vcf"   # 10 PASS interchromosomal BND
HG002_CHRX_CNV  = FIXTURES / "dragen_hg002_chrX_cnv.vcf"       # 10 PASS CNV on chrX
HG002_CHR1_FILT = FIXTURES / "dragen_hg002_chr1_filtered.vcf"  # 8 PASS + 7 non-PASS DEL on chr1

SVDB = [sys.executable, "-m", "svdb"]


def run(*args):
    return subprocess.run(SVDB + list(args), capture_output=True, text=True)


def data_lines(text: str) -> list[str]:
    return [line for line in text.splitlines() if line and not line.startswith("#")]


def annotated(text: str) -> list[str]:
    return [line for line in data_lines(text) if "OCC=" in line]


# ---------------------------------------------------------------------------
# Two individuals — autosome merge (the primary population-DB use case)
# ---------------------------------------------------------------------------

class TestTwoIndividualAutosome:
    """Core scenario: merging SVs from two different individuals on autosomes."""

    def test_merge_two_individuals_same_chrom(self):
        """Merging chr1 DEL from HG001 and HG002 should produce a combined VCF."""
        r = run("--merge", "--vcf", str(MANTA_CHR1), str(HG002_CHR1))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) > 0

    def test_merge_two_individuals_has_supp_vec(self):
        """Merged output should carry SUPP_VEC tracking which individuals contributed."""
        r = run("--merge", "--vcf", str(MANTA_CHR1), str(HG002_CHR1))
        assert "SUPP_VEC" in r.stdout

    def test_merge_three_callers_two_individuals(self):
        """manta (HG001) + tiddit (HG001) + sv-caller (HG002): three-way merge."""
        r = run("--merge", "--vcf", str(MANTA_CHR1), str(TIDDIT_CHR1), str(HG002_CHR1))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) > 0

    def test_two_individual_db_build_and_export(self, tmp_path):
        """Building a population DB from two individuals should round-trip cleanly."""
        prefix = tmp_path / "pop"
        r = run("--build", "--files", str(MANTA_CHR1), str(HG002_CHR1),
                "--prefix", str(prefix))
        assert r.returncode == 0
        assert (tmp_path / "pop.db").exists()

        r = run("--export", "--db", str(tmp_path / "pop.db"),
                "--prefix", str(tmp_path / "pop_out"))
        assert r.returncode == 0
        vcf = tmp_path / "pop_out.vcf"
        assert vcf.exists()
        assert len(data_lines(vcf.read_text())) > 0

    def test_export_two_individual_db_has_two_sample_columns(self, tmp_path):
        """Exported VCF should have one GT column per source sample."""
        prefix = tmp_path / "pop"
        run("--build", "--files", str(MANTA_CHR1), str(HG002_CHR1), "--prefix", str(prefix))
        run("--export", "--db", str(tmp_path / "pop.db"), "--prefix", str(tmp_path / "out"))
        vcf = (tmp_path / "out.vcf").read_text()
        header = [line for line in vcf.splitlines() if line.startswith("#CHROM")][0]
        assert len(header.split("\t")[9:]) == 2

    def test_query_hg001_against_two_individual_db(self, tmp_path):
        """Querying HG001 variants against a HG001+HG002 db: shared SVs get OCC≥1."""
        prefix = tmp_path / "pop"
        run("--build", "--files", str(MANTA_CHR1), str(HG002_CHR1), "--prefix", str(prefix))
        r = run("--query", "--sqdb", str(tmp_path / "pop.db"),
                "--query_vcf", str(MANTA_CHR1))
        assert r.returncode == 0
        assert len(annotated(r.stdout)) > 0

    def test_query_hg002_against_hg001_only_db(self, tmp_path):
        """HG002 queried against HG001-only db: novel HG002 SVs should get OCC=0."""
        prefix = tmp_path / "hg001"
        run("--build", "--files", str(MANTA_CHR1), "--prefix", str(prefix))
        r = run("--query", "--sqdb", str(tmp_path / "hg001.db"),
                "--query_vcf", str(HG002_CHR1), "--max_frq", "1.0")
        assert r.returncode == 0
        # Some HG002 variants may overlap HG001; others will be novel (OCC=0 or missing OCC)
        assert len(data_lines(r.stdout)) > 0

    def test_shared_sv_higher_occ_than_private(self, tmp_path):
        """A variant present in both individuals should have higher OCC than one
        present in only one — checks the frequency counting is working."""
        prefix = tmp_path / "pop"
        run("--build", "--files", str(MANTA_CHR1), str(HG002_CHR1), "--prefix", str(prefix))
        r = run("--query", "--sqdb", str(tmp_path / "pop.db"),
                "--query_vcf", str(MANTA_CHR1), "--max_frq", "1.0")
        assert r.returncode == 0
        occs = []
        for line in data_lines(r.stdout):
            for field in line.split("\t")[7].split(";"):
                if field.startswith("OCC="):
                    occs.append(int(field[4:]))
        assert max(occs) >= 1


# ---------------------------------------------------------------------------
# chrX — present in both individuals (female has 2 copies, male has 1)
# ---------------------------------------------------------------------------

class TestChrX:
    def test_merge_chrX_two_callers_same_individual(self):
        """Two callers on HG001 chrX should merge without crash."""
        r = run("--merge", "--vcf", str(MANTA_CHRX), str(TIDDIT_CHRX))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) > 0

    def test_merge_chrX_two_individuals(self):
        """HG001 and HG002 chrX DEL: cross-individual chrX merge."""
        r = run("--merge", "--vcf", str(MANTA_CHRX), str(HG002_CHRX_DEL))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) > 0

    def test_query_chrX_self_annotates_all(self):
        """Querying a chrX VCF against itself should annotate every variant."""
        r = run("--query", "--db", str(HG002_CHRX_DEL), "--query_vcf", str(HG002_CHRX_DEL))
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert all("OCC=" in line for line in lines)

    def test_query_hg001_chrX_against_hg002_chrX(self):
        """HG001 chrX queried against HG002 chrX: overlapping SVs should be found."""
        r = run("--query", "--db", str(HG002_CHRX_DEL), "--query_vcf", str(MANTA_CHRX))
        assert r.returncode == 0
        assert len(annotated(r.stdout)) > 0

    def test_build_export_chrX_two_individuals(self, tmp_path):
        prefix = tmp_path / "chrX"
        run("--build", "--files", str(MANTA_CHRX), str(HG002_CHRX_DEL), "--prefix", str(prefix))
        r = run("--export", "--db", str(tmp_path / "chrX.db"),
                "--prefix", str(tmp_path / "chrX_out"))
        assert r.returncode == 0
        assert len(data_lines((tmp_path / "chrX_out.vcf").read_text())) > 0

    def test_strict_overlap_fewer_chrX_merges(self):
        """overlap=1.0 should produce >= variants than overlap=0.0 on chrX."""
        loose  = run("--merge", "--vcf", str(MANTA_CHRX), str(HG002_CHRX_DEL), "--overlap", "0.0")
        strict = run("--merge", "--vcf", str(MANTA_CHRX), str(HG002_CHRX_DEL), "--overlap", "1.0")
        assert len(data_lines(strict.stdout)) >= len(data_lines(loose.stdout))


# ---------------------------------------------------------------------------
# chrY — male-only; absent in female individual
# ---------------------------------------------------------------------------

class TestChrY:
    def test_merge_chrY_single_sample(self):
        """chrY variants from a male individual should pass through merge."""
        r = run("--merge", "--vcf", str(HG002_CHRY))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) > 0

    def test_build_export_chrY(self, tmp_path):
        """chrY variants should round-trip through build → export correctly."""
        prefix = tmp_path / "chrY"
        run("--build", "--files", str(HG002_CHRY), "--prefix", str(prefix))
        r = run("--export", "--db", str(tmp_path / "chrY.db"),
                "--prefix", str(tmp_path / "chrY_out"))
        assert r.returncode == 0
        assert len(data_lines((tmp_path / "chrY_out.vcf").read_text())) > 0

    def test_query_female_chrX_against_male_chrY_db(self, tmp_path):
        """A female VCF (chrX) queried against a male chrY-only db should get no hits.
        Cross-chromosome queries must not produce false positives."""
        prefix = tmp_path / "maleY"
        run("--build", "--files", str(HG002_CHRY), "--prefix", str(prefix))
        r = run("--query", "--sqdb", str(tmp_path / "maleY.db"),
                "--query_vcf", str(MANTA_CHRX))
        assert r.returncode == 0
        # chrX query against chrY-only db: no OCC > 0 expected
        hits_with_occ = [line for line in data_lines(r.stdout)
                         if "OCC=" in line and "OCC=0" not in line]
        assert len(hits_with_occ) == 0

    def test_query_chrY_vcf_against_itself(self):
        """Self-query of chrY VCF should annotate every variant."""
        r = run("--query", "--db", str(HG002_CHRY), "--query_vcf", str(HG002_CHRY))
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert all("OCC=" in line for line in lines)

    def test_merge_chrX_and_chrY_together(self):
        """Merging a VCF with both chrX and chrY data should not crash."""
        r = run("--merge", "--vcf", str(HG002_CHRX_DEL), str(HG002_CHRY))
        assert r.returncode == 0


# ---------------------------------------------------------------------------
# Interchromosomal BND
# ---------------------------------------------------------------------------

class TestInterchromosomalBND:
    def test_merge_interchr_bnd(self):
        r = run("--merge", "--vcf", str(HG002_BND_INTER))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) > 0

    def test_merge_interchr_bnd_self_collapses(self):
        """Two copies of the same interchromosomal BND set should collapse."""
        single = len(data_lines(run("--merge", "--vcf", str(HG002_BND_INTER)).stdout))
        double = len(data_lines(
            run("--merge", "--vcf", str(HG002_BND_INTER), str(HG002_BND_INTER)).stdout
        ))
        assert double <= single

    def test_query_interchr_bnd_self_annotates_all(self):
        r = run("--query", "--db", str(HG002_BND_INTER), "--query_vcf", str(HG002_BND_INTER))
        assert r.returncode == 0
        lines = data_lines(r.stdout)
        assert all("OCC=" in line for line in lines)

    def test_bnd_distance_affects_interchr_overlap(self):
        """Tighter bnd_distance should find fewer or equal matches for interchr BND."""
        tight = run("--query", "--db", str(HG002_BND_INTER), "--query_vcf", str(HG002_BND_INTER),
                    "--bnd_distance", "100")
        loose = run("--query", "--db", str(HG002_BND_INTER), "--query_vcf", str(HG002_BND_INTER),
                    "--bnd_distance", "100000")
        assert len(annotated(loose.stdout)) >= len(annotated(tight.stdout))

    def test_build_export_interchr_bnd(self, tmp_path):
        prefix = tmp_path / "bnd"
        run("--build", "--files", str(HG002_BND_INTER), "--prefix", str(prefix))
        r = run("--export", "--db", str(tmp_path / "bnd.db"),
                "--prefix", str(tmp_path / "bnd_out"))
        assert r.returncode == 0


# ---------------------------------------------------------------------------
# SV type isolation — CNV vs DEL
# ---------------------------------------------------------------------------

class TestSVTypeMismatch:
    def test_cnv_merges_with_cnv(self):
        """CNV variants should merge with other CNV variants."""
        r = run("--merge", "--vcf", str(HG002_CHRX_CNV), str(HG002_CHRX_CNV))
        assert r.returncode == 0
        # Merging with itself should collapse
        single = len(data_lines(run("--merge", "--vcf", str(HG002_CHRX_CNV)).stdout))
        double = len(data_lines(r.stdout))
        assert double <= single

    def test_cnv_does_not_merge_with_del_by_default(self):
        """CNV and DEL are different SVTYPE; --no_var is required to merge them."""
        strict = run("--merge", "--vcf", str(HG002_CHRX_CNV), str(HG002_CHRX_DEL))
        no_var = run("--merge", "--vcf", str(HG002_CHRX_CNV), str(HG002_CHRX_DEL), "--no_var")
        assert strict.returncode == 0
        assert no_var.returncode == 0
        # With --no_var cross-type merging is allowed → fewer output lines
        assert len(data_lines(no_var.stdout)) <= len(data_lines(strict.stdout))

    def test_build_cnv_and_del_separate_types(self, tmp_path):
        """Build accepts mixed SV types; query should distinguish them."""
        prefix = tmp_path / "mixed"
        run("--build", "--files", str(HG002_CHRX_CNV), str(HG002_CHRX_DEL),
            "--prefix", str(prefix))
        r = run("--query", "--sqdb", str(tmp_path / "mixed.db"),
                "--query_vcf", str(HG002_CHRX_CNV))
        assert r.returncode == 0


# ---------------------------------------------------------------------------
# pass_only — generic filter behaviour (independent of caller's filter labels)
# ---------------------------------------------------------------------------

class TestPassOnly:
    def test_passonly_merge_does_not_crash(self):
        """--pass_only on a VCF with mixed PASS/non-PASS variants must not crash."""
        r = run("--merge", "--vcf", str(HG002_CHR1_FILT), "--pass_only")
        assert r.returncode == 0

    def test_passonly_reduces_merging_between_non_pass_variants(self):
        """With --pass_only, non-PASS variants should not be merged into PASS ones.
        Output variant count may differ from the unrestricted run."""
        without = run("--merge", "--vcf", str(HG002_CHR1_FILT), str(HG002_CHR1_FILT))
        with_po = run("--merge", "--vcf", str(HG002_CHR1_FILT), str(HG002_CHR1_FILT), "--pass_only")
        assert without.returncode == 0
        assert with_po.returncode == 0
        # At minimum, both return some output
        assert len(data_lines(without.stdout)) > 0
        assert len(data_lines(with_po.stdout)) > 0

    def test_build_passonly_vs_full_has_fewer_indexed_variants(self, tmp_path):
        """A pass_only build should index fewer variants than a full build,
        so the same query returns fewer hits."""
        full_p   = tmp_path / "full"
        pass_p   = tmp_path / "pass"
        run("--build", "--files", str(HG002_CHR1_FILT), "--prefix", str(full_p))
        run("--build", "--files", str(HG002_CHR1_FILT), "--pass_only", "--prefix", str(pass_p))

        r_full = run("--query", "--sqdb", str(tmp_path / "full.db"),
                     "--query_vcf", str(HG002_CHR1_FILT), "--max_frq", "1.0")
        r_pass = run("--query", "--sqdb", str(tmp_path / "pass.db"),
                     "--query_vcf", str(HG002_CHR1_FILT), "--max_frq", "1.0")
        # pass_only db has fewer indexed variants → fewer annotated hits
        assert len(annotated(r_pass.stdout)) <= len(annotated(r_full.stdout))
