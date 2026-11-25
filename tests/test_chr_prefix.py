"""Issue #115: chromosome names with and without a "chr" prefix must match.

Fixtures (tests/fixtures/chr_prefix): the same four variants (two DELs, a 1->2 BND and a
mitochondrial DEL) written once as 1/2/MT (sample S1) and once as chr1/chr2/chrM (sample S2).
"""

import subprocess
import sys
from pathlib import Path

import pytest

from svdb.read_vcf import readVCFLine
from svdb.vcf_utils import chrom_key

FIXTURES = Path(__file__).parent / "fixtures" / "chr_prefix"
NOCHR = FIXTURES / "nochr.vcf"
CHR = FIXTURES / "chr.vcf"
IDS = {"del1", "del2", "bnd1", "mt1"}

SVDB = [sys.executable, "-m", "svdb"]


def run(*args):
    r = subprocess.run(SVDB + [str(a) for a in args], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def records(vcf_text):
    out = []
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        out.append({"CHROM": f[0], "ID": f[2], "INFO": info})
    return out


def build(vcfs, prefix):
    run("--build", "--files", *vcfs, "--prefix", prefix)
    return prefix.with_suffix(".db")


@pytest.mark.parametrize("name, key", [
    ("1", "1"), ("chr1", "1"), ("Chr1", "1"), ("CHR1", "1"),
    ("X", "X"), ("chrX", "X"),
    ("M", "MT"), ("MT", "MT"), ("chrM", "MT"), ("chrMT", "MT"),
    ("chrUn_KI270302v1", "Un_KI270302v1"),
    ("HLA-A*01:01", "HLA-A*01:01"),
    ("scaffold_chr1", "scaffold_chr1"),
])
def test_chrom_key(name, key):
    assert chrom_key(name) == key


def test_parse_keeps_names_and_sets_keys():
    v = readVCFLine("chr1\t3000000\tbnd1\tN\tN]2:5000000]\t.\tPASS\tSVTYPE=BND\tGT\t0/1")
    assert (v.chrA, v.chrB) == ("chr1", "2")
    assert (v.chrA_key, v.chrB_key) == ("1", "2")


@pytest.mark.parametrize("first, second", [(NOCHR, CHR), (CHR, NOCHR)])
def test_merge_matches_across_styles(first, second):
    merged = records(run("--merge", "--vcf", first, second))
    assert len(merged) == len(IDS)
    assert all(r["INFO"]["set"] == "Intersection" for r in merged)


def test_merge_keeps_name_of_kept_record():
    merged = records(run("--merge", "--vcf", f"{CHR}:chr", f"{NOCHR}:nochr", "--priority", "chr,nochr"))
    assert {r["CHROM"] for r in merged} == {"chr1", "chrM"}


@pytest.mark.parametrize("query, db", [(CHR, NOCHR), (NOCHR, CHR)])
def test_query_db_matches_across_styles(query, db):
    hits = records(run("--query", "--query_vcf", query, "--db", db))
    assert {r["ID"]: r["INFO"].get("OCC") for r in hits} == dict.fromkeys(IDS, "1")


def test_query_bedpedb_matches_across_styles():
    hits = records(run("--query", "--query_vcf", NOCHR, "--bedpedb", FIXTURES / "chr.bedpe"))
    assert {r["ID"]: r["INFO"].get("OCC") for r in hits} == dict.fromkeys(IDS, "3")


@pytest.mark.parametrize("query, db_vcf", [(CHR, NOCHR), (NOCHR, CHR)])
def test_query_sqdb_matches_across_styles(query, db_vcf, tmp_path):
    db = build([db_vcf], tmp_path / "svdb")
    hits = records(run("--query", "--query_vcf", query, "--sqdb", db))
    assert {r["ID"]: r["INFO"].get("OCC") for r in hits} == dict.fromkeys(IDS, "1")


def test_export_clusters_across_styles(tmp_path):
    db = build([NOCHR, CHR], tmp_path / "svdb")
    run("--export", "--db", db, "--prefix", tmp_path / "out")
    exported = records((tmp_path / "out.vcf").read_text())
    assert len(exported) == len(IDS)
    assert all(r["INFO"]["OCC"] == "2" for r in exported)
