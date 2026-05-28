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
INS_SVLEN_RANGE_DIR = FIXTURES / "ins_svlen_range"
MANTA = FIXTURES / "manta_chr1_del.vcf"
TIDDIT = FIXTURES / "tiddit_chr1_del.vcf"
TRUTH = FIXTURES / "truth_chr1_del.vcf"
MANTA_BND = FIXTURES / "manta_bnd.vcf"
CNVKIT = FIXTURES / "cnvkit_chr1_del.vcf"

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

    def test_build_pass_only(self, tmp_path):
        prefix = tmp_path / "svdb_pass"
        r = run("--build", "--files", str(MANTA), "--pass_only", "--prefix", str(prefix))
        assert r.returncode == 0
        assert (tmp_path / "svdb_pass.db").exists()

    def test_build_passonly_deprecated_alias(self, tmp_path):
        prefix = tmp_path / "svdb_pass_alias"
        r = run("--build", "--files", str(MANTA), "--passonly", "--prefix", str(prefix))
        assert r.returncode == 0
        assert (tmp_path / "svdb_pass_alias.db").exists()
        assert "deprecated" in r.stderr.lower()

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

    def test_export_coarse(self, db, tmp_path):
        prefix = tmp_path / "coarse"
        r = run("--export", "--db", str(db), "--coarse", "--epsilon", "500",
                "--min_pts", "2", "--prefix", str(prefix))
        assert r.returncode == 0

    def test_export_dbscan_deprecated_alias(self, db, tmp_path):
        prefix = tmp_path / "dbscan_alias"
        r = run("--export", "--db", str(db), "--DBSCAN", "--prefix", str(prefix))
        assert r.returncode == 0
        assert "deprecated" in r.stderr.lower()

    def test_export_memory_flag(self, db, tmp_path):
        prefix = tmp_path / "mem"
        r = run("--export", "--db", str(db), "--memory", "--prefix", str(prefix))
        assert r.returncode == 0


# ---------------------------------------------------------------------------
# Export — strip_chr
# ---------------------------------------------------------------------------

_CHR_VCF_CONTENT = """\
##fileformat=VCFv4.1
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t1000\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000\tGT\t0/1
"""


class TestExportStripChr:
    @pytest.fixture
    def chr_db(self, tmp_path):
        vcf = tmp_path / "chr_sample.vcf"
        vcf.write_text(_CHR_VCF_CONTENT)
        prefix = tmp_path / "svdb"
        run("--build", "--files", str(vcf), "--prefix", str(prefix))
        return tmp_path / "svdb.db"

    def test_chr_prefix_preserved_without_flag(self, chr_db, tmp_path):
        """Without --strip_chr, exported CHROM column keeps the chr prefix."""
        prefix = tmp_path / "out"
        r = run("--export", "--db", str(chr_db), "--prefix", str(prefix))
        assert r.returncode == 0
        lines = vcf_data_lines((tmp_path / "out.vcf").read_text())
        assert all(line.startswith("chr") for line in lines)

    def test_strip_chr_removes_prefix(self, chr_db, tmp_path):
        """With --strip_chr, exported CHROM column has no chr prefix."""
        prefix = tmp_path / "out"
        r = run("--export", "--db", str(chr_db), "--strip_chr", "--prefix", str(prefix))
        assert r.returncode == 0
        lines = vcf_data_lines((tmp_path / "out.vcf").read_text())
        assert len(lines) > 0
        assert all(not line.startswith("chr") for line in lines)


# ---------------------------------------------------------------------------
# Query: exported VCF round-trip for insertions with SVLEN variation
# ---------------------------------------------------------------------------

def _ins_vcf(sample_name: str, pos: int, seq: str) -> str:
    """Minimal single-sample INS VCF with an actual sequence in ALT."""
    svlen = len(seq)
    return (
        "##fileformat=VCFv4.1\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        f"chr1\t{pos}\tins1\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={svlen}\tGT\t0/1\n"
    )


class TestQueryInsVcfRoundTrip:
    """A sample that contributed to a DB cluster must be findable when
    querying the exported VCF, even if its SVLEN differs from the cluster
    representative's SVLEN.

    Scenario: 3 samples, all INS at chr1:1000.
      sample_a: SVLEN=64  (8-copy tandem repeat; appears twice → becomes representative)
      sample_b: SVLEN=64
      sample_c: SVLEN=80  (10-copy; ratio vs representative = 64/80 = 0.80 < default 0.90)

    After build → export (permissive ratio=0.79) → query sample_c with default
    ratio=0.90 against the exported VCF: the cluster representative is SVLEN=64
    but sample_c (SVLEN=80) is in the cluster.  Without the per-member SVLEN fix
    the query returns no OCC for sample_c.
    """

    # tandem-repeat sequences: 8-copy and 10-copy of the same unit
    _UNIT = "ATCGATCG"
    _SEQ_64 = _UNIT * 8                             # 64 bp  (representative)
    _SEQ_80 = _UNIT * 10                            # 80 bp  (outlier, ratio 64/80=0.80 < 0.90)

    @pytest.fixture
    def round_trip(self, tmp_path):
        """Build DB from 3 samples, export to VCF, return (exported_vcf, query_vcf_c)."""
        vcf_a = tmp_path / "sample_a.vcf"
        vcf_b = tmp_path / "sample_b.vcf"
        vcf_c = tmp_path / "sample_c.vcf"
        vcf_a.write_text(_ins_vcf("sample_a", 1000, self._SEQ_64))
        vcf_b.write_text(_ins_vcf("sample_b", 1000, self._SEQ_64))
        vcf_c.write_text(_ins_vcf("sample_c", 1000, self._SEQ_80))

        prefix = tmp_path / "test_db"
        r = run(
            "--build", "--files", str(vcf_a), str(vcf_b), str(vcf_c),
            "--prefix", str(prefix),
        )
        assert r.returncode == 0, r.stderr

        export_prefix = tmp_path / "exported"
        r = run(
            "--export", "--db", str(tmp_path / "test_db.db"),
            "--prefix", str(export_prefix),
            # 64/80=0.80 ≥ 0.79 → all three land in one cluster;
            # sequence similarity is ~0.80 (≥ default 0.75) so no --ins_seq_similarity override needed
            "--ins_svlen_ratio", "0.79",
        )
        assert r.returncode == 0, r.stderr

        return tmp_path / "exported.vcf", vcf_c

    def test_svlen_outlier_member_annotated_by_range_check(self, round_trip, tmp_path):
        """Querying sample_c (SVLEN=80) against a cluster whose representative
        is SVLEN=64 must return OCC > 0 (default ins_svlen_ratio=0.90)."""
        exported_vcf, query_vcf = round_trip
        result = run(
            "--query",
            "--query_vcf", str(query_vcf),
            "--db", str(exported_vcf),
        )
        assert result.returncode == 0, result.stderr

        data_lines = vcf_data_lines(result.stdout)
        assert data_lines, "no variant lines in query output"

        # Every INS line must have OCC annotated (not missing)
        unannotated = [
            ln for ln in data_lines
            if "INS" in ln and "OCC=" not in ln
        ]
        assert not unannotated, (
            f"{len(unannotated)} INS variant(s) missing OCC — "
            "cluster member not matched against its own cluster in exported VCF"
        )


class TestQueryInsVcfRoundTripFixture:
    """Real-world regression test using 9 Sniffles2 VCFs at chrX:128937049.

    All 9 samples share an INS at this locus with SVLENs ranging from 496–602 bp.
    After build → export, the cluster representative has SVLEN=500 (HG01312).
    Querying HG04214 (SVLEN=564) with --db against the exported VCF previously
    returned no OCC because 500/564=0.887 < default ins_svlen_ratio=0.90.
    With per-member SVLEN stored in VARIANTS, the range check passes.
    """

    @pytest.fixture
    def tiny_round_trip(self, tmp_path):
        vcfs = sorted(INS_SVLEN_RANGE_DIR.glob("*.vcf"))
        assert vcfs, f"tiny fixture directory empty: {INS_SVLEN_RANGE_DIR}"

        r = run("--build", "--files", *[str(v) for v in vcfs], "--prefix", str(tmp_path / "db"))
        assert r.returncode == 0, r.stderr

        r = run("--export", "--db", str(tmp_path / "db.db"), "--prefix", str(tmp_path / "exp"))
        assert r.returncode == 0, r.stderr

        return tmp_path / "exp.vcf"

    def test_hg04214_svlen_564_matched_against_cluster_rep_500(self, tiny_round_trip):
        """HG04214 (SVLEN=564) must be annotated against a cluster whose representative
        is SVLEN=500 (500/564=0.887 < default 0.90 without the range fix)."""
        query_vcf = INS_SVLEN_RANGE_DIR / "HG04214.vcf"
        result = run("--query", "--query_vcf", str(query_vcf), "--db", str(tiny_round_trip))
        assert result.returncode == 0, result.stderr

        data_lines = vcf_data_lines(result.stdout)
        assert data_lines, "no variant lines in query output"

        unannotated = [ln for ln in data_lines if "INS" in ln and "OCC=" not in ln]
        assert not unannotated, (
            f"{len(unannotated)} INS variant(s) missing OCC — "
            "SVLEN-range fix not working for tiny fixture"
        )

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


# ---------------------------------------------------------------------------
# CNVkit parenthesised CIPOS/CIEND preservation (issue #72)
# ---------------------------------------------------------------------------


class TestCNVkitCiposPreservation:
    """Regression tests for issue #72.

    CNVkit emits confidence intervals as CIPOS=(0,166417) — parenthesised,
    with a zero lower bound.  The bug: SVDB was converting these to
    CIPOS=.,166417 (lower bound replaced with the VCF missing-value sentinel).

    Tests run the full merge and build→export pipelines and assert that no
    CIPOS or CIEND value in the output contains '.' as a bound.  This catches
    the corruption wherever it occurs in the pipeline without being tied to
    a specific code location.
    """

    @staticmethod
    def _ci_values(vcf_text: str) -> list[str]:
        """Return every raw CIPOS=... and CIEND=... token from data lines."""
        values = []
        for line in vcf_text.splitlines():
            if not line or line.startswith("#"):
                continue
            for part in line.split("\t")[7].split(";"):
                if part.startswith("CIPOS=") or part.startswith("CIEND="):
                    values.append(part)
        return values

    def test_merge_preserves_cipos_lower_bound(self):
        """After merge the CIPOS lower bound must be a number, not '.'."""
        r = run("--merge", "--vcf", str(CNVKIT))
        assert r.returncode == 0
        ci_tokens = self._ci_values(r.stdout)
        assert ci_tokens, "expected CIPOS/CIEND tokens in merged output"
        bad = [t for t in ci_tokens if t.split("=", 1)[1].startswith(".")]
        assert not bad, (
            "CIPOS/CIEND lower bound corrupted to '.' in merge output:\n"
            + "\n".join(bad)
        )

    def test_build_export_preserves_cipos_bounds(self, tmp_path):
        """After build→export neither CIPOS bound must be '.'."""
        prefix = tmp_path / "svdb"
        r_build = run("--build", "--files", str(CNVKIT), "--prefix", str(prefix))
        assert r_build.returncode == 0

        out_prefix = tmp_path / "export"
        r_export = run("--export", "--db", str(tmp_path / "svdb.db"),
                       "--prefix", str(out_prefix))
        assert r_export.returncode == 0

        exported_vcf = (tmp_path / "export.vcf").read_text()
        ci_tokens = self._ci_values(exported_vcf)
        assert ci_tokens, "expected CIPOS/CIEND tokens in exported output"
        bad = [t for t in ci_tokens if "." in t.split("=", 1)[1].split(",")]
        assert not bad, (
            "CIPOS/CIEND bound corrupted to '.' in build→export output:\n"
            + "\n".join(bad)
        )


# ---------------------------------------------------------------------------
# Query consistency — --db vs --sqdb (issue #70)
# ---------------------------------------------------------------------------


class TestQueryDbVsSqDbConsistency:
    """Regression tests for issue #70: --db and --sqdb must produce identical
    OCC/FRQ annotations when the SQLite database was built from the same VCF
    used as the direct --db input.

    These tests document the expected behaviour.  If they fail, the two code
    paths have drifted and the discrepancy needs investigation before fixing.
    """

    @pytest.fixture
    def sqdb_from_manta(self, tmp_path):
        """SQLite DB built from MANTA only — mirrors using MANTA as --db."""
        prefix = tmp_path / "svdb"
        run("--build", "--files", str(MANTA), "--prefix", str(prefix))
        return tmp_path / "svdb.db"

    @staticmethod
    def _annotations(vcf_text: str, occ_tag: str = "OCC", frq_tag: str = "FRQ") -> dict:
        """Return {(chrom, pos): (occ, frq)} for every annotated variant."""
        result = {}
        for line in vcf_text.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, pos, info = fields[0], fields[1], fields[7]
            occ = frq = None
            for part in info.split(";"):
                if part.startswith(occ_tag + "="):
                    occ = int(part.split("=", 1)[1])
                elif part.startswith(frq_tag + "="):
                    frq = float(part.split("=", 1)[1])
            if occ is not None:
                result[(chrom, pos)] = (occ, frq)
        return result

    def test_same_variants_annotated(self, sqdb_from_manta):
        """Every variant annotated by --db should also be annotated by --sqdb
        built from the same source, and vice versa."""
        vcf_r = run("--query", "--db", str(MANTA), "--query_vcf", str(TRUTH))
        sq_r  = run("--query", "--sqdb", str(sqdb_from_manta), "--query_vcf", str(TRUTH))
        assert vcf_r.returncode == 0
        assert sq_r.returncode == 0

        vcf_hits  = set(self._annotations(vcf_r.stdout).keys())
        sqdb_hits = set(self._annotations(sq_r.stdout).keys())

        only_vcf  = vcf_hits - sqdb_hits
        only_sqdb = sqdb_hits - vcf_hits
        assert vcf_hits == sqdb_hits, (
            f"Annotated variant sets differ.\n"
            f"  Only in --db:   {sorted(only_vcf)}\n"
            f"  Only in --sqdb: {sorted(only_sqdb)}"
        )

    def test_occ_values_match(self, sqdb_from_manta):
        """OCC counts must be identical for every variant annotated by both paths."""
        vcf_r = run("--query", "--db", str(MANTA), "--query_vcf", str(TRUTH))
        sq_r  = run("--query", "--sqdb", str(sqdb_from_manta), "--query_vcf", str(TRUTH))

        vcf_ann  = self._annotations(vcf_r.stdout)
        sqdb_ann = self._annotations(sq_r.stdout)

        mismatches = [
            f"  {key}: --db OCC={vcf_ann[key][0]}, --sqdb OCC={sqdb_ann[key][0]}"
            for key in vcf_ann.keys() & sqdb_ann.keys()
            if vcf_ann[key][0] != sqdb_ann[key][0]
        ]
        assert not mismatches, (
            "OCC mismatch between --db and --sqdb:\n" + "\n".join(mismatches)
        )

    def test_frq_values_match(self, sqdb_from_manta):
        """FRQ values must be identical for every variant annotated by both paths."""
        vcf_r = run("--query", "--db", str(MANTA), "--query_vcf", str(TRUTH))
        sq_r  = run("--query", "--sqdb", str(sqdb_from_manta), "--query_vcf", str(TRUTH))

        vcf_ann  = self._annotations(vcf_r.stdout)
        sqdb_ann = self._annotations(sq_r.stdout)

        mismatches = [
            f"  {key}: --db FRQ={vcf_ann[key][1]}, --sqdb FRQ={sqdb_ann[key][1]}"
            for key in vcf_ann.keys() & sqdb_ann.keys()
            if vcf_ann[key][1] != sqdb_ann[key][1]
        ]
        assert not mismatches, (
            "FRQ mismatch between --db and --sqdb:\n" + "\n".join(mismatches)
        )
