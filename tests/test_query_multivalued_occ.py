"""Integration tests for issue #105: --query with --in_occ/--in_frq must not crash on
database records whose OCC/FRQ tag holds more than one value, and such records must
not be reported as hits.

gnomAD-SV v2.1 multi-allelic CNVs (SVTYPE=MCNV) carry per-copy-state lists, e.g.
AC=1731,6358,2402,... These are counts of individuals per copy state (CN=2 included),
not carrier counts, so there is no single OCC/FRQ value to report for them.

Fixtures (tests/fixtures/mcnv_no_var), one case per chromosome:
  chr1  DEL overlapping a real gnomAD-SV MCNV record (the reproducer from the issue)
  chr2  DEL with two db candidates: a biallelic DEL (AC=11) and an MCNV that overlaps
        the query better
  chr3  INS whose distance to both ends of a 200 bp MCNV is within --ins_distance
"""

import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "mcnv_no_var"
QUERY = FIXTURES / "query.vcf"
DB = FIXTURES / "db.vcf"

SVDB = [sys.executable, "-m", "svdb"]
QUERY_OPTS = ["--query", "--query_vcf", str(QUERY), "--db", str(DB),
              "--in_occ", "AC", "--in_frq", "AF", "--ins_distance", "500"]


def run(*args):
    return subprocess.run(SVDB + list(args), capture_output=True, text=True)


def info_by_id(vcf_text):
    """Map record ID -> parsed INFO dict."""
    out = {}
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        out[f[2]] = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
    return out


@pytest.fixture(scope="module")
def no_var():
    return run(*QUERY_OPTS, "--no_var")


@pytest.fixture(scope="module")
def default():
    return run(*QUERY_OPTS)


class TestNoVarSkipsMultiValuedRecords:

    def test_query_succeeds(self, no_var):
        assert no_var.returncode == 0, no_var.stderr

    def test_all_query_records_written(self, no_var):
        assert no_var.returncode == 0, no_var.stderr
        assert set(info_by_id(no_var.stdout)) == {
            "del_vs_mcnv", "del_vs_del_and_mcnv", "ins_vs_short_mcnv"}

    def test_del_overlapping_mcnv_is_not_annotated(self, no_var):
        assert no_var.returncode == 0, no_var.stderr
        assert "OCC" not in info_by_id(no_var.stdout)["del_vs_mcnv"]

    def test_biallelic_del_wins_over_better_overlapping_mcnv(self, no_var):
        assert no_var.returncode == 0, no_var.stderr
        info = info_by_id(no_var.stdout)["del_vs_del_and_mcnv"]
        assert int(info["OCC"]) == 11
        assert float(info["FRQ"]) == pytest.approx(0.011)

    def test_ins_near_short_mcnv_is_not_annotated(self, no_var):
        assert no_var.returncode == 0, no_var.stderr
        assert "OCC" not in info_by_id(no_var.stdout)["ins_vs_short_mcnv"]


class TestDefaultModeControl:
    """Without --no_var no query matches an MCNV record; this already works."""

    def test_query_succeeds(self, default):
        assert default.returncode == 0, default.stderr

    def test_annotations(self, default):
        info = info_by_id(default.stdout)
        assert "OCC" not in info["del_vs_mcnv"]
        assert int(info["del_vs_del_and_mcnv"]["OCC"]) == 11
        assert float(info["del_vs_del_and_mcnv"]["FRQ"]) == pytest.approx(0.011)
        assert "OCC" not in info["ins_vs_short_mcnv"]
