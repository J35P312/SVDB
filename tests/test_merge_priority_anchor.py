"""Integration tests for issue #102: the merged line must be the highest-priority
record of its cluster, so ALT and the priority-only INFO tags (SVTYPE, END, SVLEN)
come from the same record.

Fixtures (tests/fixtures/tra_bnd_priority) are the reproducer from the issue:
a Sawfish-style breakend pair (bracket ALT, no END) and a DeBreak-style <TRA>
(CHR2 + END) 200 bp to its left, plus a DEL pair where the two callers are 10 bp apart.
"""

import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "tra_bnd_priority"
SAWFISH = FIXTURES / "sawfish.vcf"
DEBREAK = FIXTURES / "debreak.vcf"

SVDB = [sys.executable, "-m", "svdb"]
MERGE_OPTS = ["--overlap", "0.7", "--bnd_distance", "2000", "--ins_distance", "100"]


def run(*args):
    return subprocess.run(SVDB + list(args), capture_output=True, text=True)


def records(vcf_text):
    """Return data lines as dicts with CHROM, POS, ID, ALT and parsed INFO."""
    out = []
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
        out.append({"CHROM": f[0], "POS": int(f[1]), "ID": f[2], "ALT": f[4], "INFO": info})
    return out


def merge(*vcfs, priority=None):
    args = ["--merge", "--vcf", *vcfs]
    if priority:
        args += ["--priority", priority]
    r = run(*args, *MERGE_OPTS)
    assert r.returncode == 0, r.stderr
    return r.stdout


def find(recs, id_prefix):
    hits = [r for r in recs if r["ID"].startswith(id_prefix)]
    assert len(hits) == 1, [r["ID"] for r in recs]
    return hits[0]


@pytest.fixture
def sawfish_first():
    return merge(f"{SAWFISH}:sawfish", f"{DEBREAK}:debreak", priority="sawfish,debreak")


@pytest.fixture
def debreak_first():
    return merge(f"{SAWFISH}:sawfish", f"{DEBREAK}:debreak", priority="debreak,sawfish")


class TestAnchorIsHighestPriority:

    def test_bnd_cluster_written_from_sawfish_breakend(self, sawfish_first):
        rec = find(records(sawfish_first), "sawfish_bnd")
        assert rec["CHROM"] == "chr4"
        assert rec["POS"] == 3000200
        assert "TRA" not in rec["ALT"]
        assert rec["INFO"]["SVTYPE"] == "BND"
        assert rec["INFO"]["set"] == "Intersection"

    def test_bnd_cluster_written_from_debreak_tra(self, debreak_first):
        rec = find(records(debreak_first), "debreak_tra:")
        assert rec["POS"] == 3000000
        assert rec["ALT"] == "<TRA>"
        assert rec["INFO"]["END"] == "5000005"
        assert rec["INFO"]["set"] == "Intersection"

    def test_del_cluster_written_from_sawfish(self, sawfish_first):
        rec = find(records(sawfish_first), "sawfish_del")
        assert rec["POS"] == 1000000
        assert rec["INFO"]["END"] == "1005000"

    def test_del_cluster_written_from_debreak_when_it_is_right_of_sawfish(self, debreak_first):
        # debreak_del sits 10 bp right of sawfish_del; with debreak ranked first it
        # must still be the written line, so the match has to be found to its left.
        rec = find(records(debreak_first), "debreak_del")
        assert rec["POS"] == 1000010
        assert rec["INFO"]["END"] == "1005010"
        assert rec["INFO"]["set"] == "Intersection"

    def test_without_priority_file_order_decides_anchor(self):
        recs = records(merge(str(SAWFISH), str(DEBREAK)))
        rec = find(recs, "sawfish_bnd")
        assert rec["POS"] == 3000200
        assert "TRA" not in rec["ALT"]

    @pytest.mark.parametrize("order", ["sawfish,debreak", "debreak,sawfish"])
    def test_cluster_membership_independent_of_priority(self, order):
        recs = records(merge(f"{SAWFISH}:sawfish", f"{DEBREAK}:debreak", priority=order))
        assert len(recs) == 3
        assert sorted(r["INFO"]["set"] for r in recs) == ["Intersection", "Intersection", "debreak"]


class TestMergedOutputIsQueryable:

    @pytest.mark.parametrize("order", ["sawfish,debreak", "debreak,sawfish"])
    def test_every_tra_record_has_end(self, order):
        recs = records(merge(f"{SAWFISH}:sawfish", f"{DEBREAK}:debreak", priority=order))
        for rec in recs:
            if "TRA" in rec["ALT"]:
                assert "END" in rec["INFO"], rec["ID"]

    def test_query_merged_against_sawfish(self, sawfish_first, tmp_path):
        merged = tmp_path / "merged.vcf"
        merged.write_text(sawfish_first)
        r = run("--query", "--query_vcf", str(merged), "--db", str(SAWFISH),
                "--prefix", str(tmp_path / "q"))
        assert r.returncode == 0, r.stderr
