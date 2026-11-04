"""Issue #111: --svtype_match {exact,hierarchical} for merge and query.

hierarchical: a type also matches its parent or child (VCF 4.4 colon subtypes;
CNV covers DEL and DUP), never its siblings.

Fixtures (tests/fixtures/svtype_match): real gnomAD-SV v4.1 chr1 records, and
DRAGEN, Delly cnv, Sawfish and MELT style records at the same loci.
"""

import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "svtype_match"
GNOMAD = FIXTURES / "gnomad_v4.1_sites.vcf"
DRAGEN = FIXTURES / "dragen.vcf"
DELLY = FIXTURES / "delly_cnv.vcf"
SAWFISH = FIXTURES / "sawfish.vcf"
MELT = FIXTURES / "melt.vcf"

SVDB = [sys.executable, "-m", "svdb"]
MODES = ["exact", "hierarchical"]


def run(*args, check=True):
    r = subprocess.run(SVDB + [str(a) for a in args], capture_output=True, text=True)
    if check:
        assert r.returncode == 0, r.stderr
    return r


def records(vcf_text):
    out = []
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
        out.append({"POS": int(f[1]), "ID": f[2], "INFO": info})
    return out


def occ_at(vcf_text, pos):
    hits = [r for r in records(vcf_text) if r["POS"] == pos]
    assert len(hits) == 1
    return float(hits[0]["INFO"].get("OCC", 0))


def merged_clusters(vcf_text):
    """Source IDs of every merged record (merged IDs are 'id:file|id:file')."""
    return {frozenset(x.rsplit(":", 1)[0] for x in r["ID"].split("|"))
            for r in records(vcf_text) if "|" in r["ID"]}


@pytest.mark.parametrize("a, b, exact, hierarchical", [
    ("DEL", "DEL", True, True),
    ("INS:ME:ALU", "INS", False, True),
    ("INS", "INS:ME:ALU", False, True),
    ("INS:ME:ALU", "INS:ME", False, True),
    ("INS:ME:ALU", "INS:ME:LINE1", False, False),
    ("CNV", "DEL", False, True),
    ("DUP", "CNV", False, True),
    ("DEL:ME:LINE1", "CNV", False, True),
    ("CNV:TR", "CNV", False, True),
    ("CNV:TR", "DEL", False, False),
    ("DEL", "DUP", False, False),
    ("INV", "DEL", False, False),
    ("BND", "BND", True, True),
])
def test_svtypes_match(a, b, exact, hierarchical):
    from svdb.svtype import svtypes_match
    assert svtypes_match(a, b, hierarchical=False) is exact
    assert svtypes_match(a, b, hierarchical=True) is hierarchical


# (query vcf, position, hit with exact, hit with hierarchical) against gnomAD
CASES = [
    pytest.param(DRAGEN, 46000, False, True, id="DEL~CNV"),
    pytest.param(DRAGEN, 11000, False, False, id="DEL!~DUP"),
    pytest.param(DRAGEN, 79086, False, True, id="DEL~DEL:ME:LINE1"),
    pytest.param(DELLY, 11000, False, True, id="CNV~DUP"),
    pytest.param(SAWFISH, 61300, False, True, id="INS~INS:ME:ALU"),
    pytest.param(MELT, 61305, True, True, id="ALU~ALU"),
    pytest.param(MELT, 634600, False, False, id="ALU!~LINE1"),
]


@pytest.fixture(scope="module")
def gnomad_sqdb(tmp_path_factory):
    prefix = tmp_path_factory.mktemp("svtype_match") / "gnomad"
    run("--build", "--files", GNOMAD, "--prefix", prefix)
    return Path(f"{prefix}.db")


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("query, pos, exact, hierarchical", CASES)
def test_query_sqdb(gnomad_sqdb, query, pos, exact, hierarchical, mode):
    out = run("--query", "--sqdb", gnomad_sqdb, "--query_vcf", query, "--svtype_match", mode).stdout
    expected = exact if mode == "exact" else hierarchical
    assert (occ_at(out, pos) > 0) is expected


# gnomAD CNV records carry no AC, so --db skips them: DEL~CNV is tested against DRAGEN below
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("query, pos, exact, hierarchical", [c for c in CASES if c.id != "DEL~CNV"])
def test_query_db_gnomad(query, pos, exact, hierarchical, mode):
    out = run("--query", "--db", GNOMAD, "--in_occ", "AC", "--in_frq", "AF",
              "--query_vcf", query, "--svtype_match", mode).stdout
    expected = exact if mode == "exact" else hierarchical
    assert (occ_at(out, pos) > 0) is expected


@pytest.mark.parametrize("mode, expected", [("exact", False), ("hierarchical", True)])
def test_query_db_cnv_hits_del(mode, expected):
    out = run("--query", "--db", DRAGEN, "--query_vcf", DELLY, "--svtype_match", mode).stdout
    assert (occ_at(out, 11000) > 0) is expected


@pytest.mark.parametrize("mode, expected", [
    ("exact", set()),
    ("hierarchical", {
        frozenset({"DRAGEN:LOSS:chr1:46001-99350", "gnomAD-SV_v3_CNV_chr1_e3665ac3"}),
        frozenset({"DRAGEN:DEL:12011:0:1:0:0:0", "gnomAD-SV_v3_DEL_chr1_733c4ef0"}),
    }),
])
def test_merge_del_cnv(mode, expected):
    out = run("--merge", "--vcf", DRAGEN, GNOMAD, "--svtype_match", mode).stdout
    assert merged_clusters(out) == expected


@pytest.mark.parametrize("mode, expected", [
    ("exact", {frozenset({"ALU_umary_1", "gnomAD-SV_v3_INS_chr1_ce3bb71c"})}),
    ("hierarchical", {frozenset({"sawfish:0:118:0:0", "ALU_umary_1", "gnomAD-SV_v3_INS_chr1_ce3bb71c"})}),
])
def test_merge_insertions(mode, expected):
    out = run("--merge", "--vcf", SAWFISH, MELT, GNOMAD, "--svtype_match", mode).stdout
    assert merged_clusters(out) == expected


@pytest.mark.parametrize("module", [
    ["--query", "--db", GNOMAD, "--query_vcf", DRAGEN],
    ["--merge", "--vcf", DRAGEN, GNOMAD],
])
def test_no_var_conflicts_with_hierarchical(module):
    r = run(*module, "--no_var", "--svtype_match", "hierarchical", check=False)
    assert r.returncode != 0
    assert "not allowed with" in r.stderr
