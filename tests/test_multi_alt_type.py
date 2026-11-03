"""Issue #107: SV type of multi-ALT symbolic records.

Fixtures (tests/fixtures/multi_alt_type) follow the VCF layout of the callers that
emit such records; each header names the caller and version it is modelled on.
Loci on GRCh37 chr1: 46000-111250 (MCNV), 16800001-16860000 (DUP),
20100001-20150000 (DEL), 152555001-152590000 (mixed / LOH).
"""

import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from svdb.read_vcf import readVCFLine

FIXTURES = Path(__file__).parent / "fixtures" / "multi_alt_type"
GCNV_SEGMENTS = FIXTURES / "gcnv_segments_v4.1.4.vcf"
GCNV_COHORT = FIXTURES / "gcnv_joint_cohort.vcf"
DRAGEN = FIXTURES / "dragen_cnv.vcf"
DRAGEN_ASCN = FIXTURES / "dragen_ascn_legacy.vcf"
CANVAS = FIXTURES / "canvas_cohort.vcf"
GNOMAD = FIXTURES / "gnomad_v2.1_mcnv.vcf"

SVDB = [sys.executable, "-m", "svdb"]


def run(*args):
    r = subprocess.run(SVDB + [str(a) for a in args], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def parse(vcf, record_id):
    for line in vcf.read_text().splitlines():
        if not line.startswith("#") and line.split("\t")[2] == record_id:
            return readVCFLine(line)
    raise KeyError(record_id)


def records(vcf_text):
    out = []
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
        out.append({"POS": int(f[1]), "ID": f[2], "ALT": f[4], "INFO": info})
    return out


@pytest.mark.parametrize("vcf, record_id, expected", [
    # gnomAD-SV v2.1 MCNV, sites only
    (GNOMAD, "gnomAD-SV_v2.1_MCNV_1_1", "CNV"),
    (GNOMAD, "gnomAD-SV_v2.1_MCNV_1_2", "CNV"),
    # gCNV <= 4.1.4: every record <DEL>,<DUP>, haploid GT picks the allele
    (GCNV_SEGMENTS, "CNV_1_20100001_20150000", "DEL"),
    (GCNV_SEGMENTS, "CNV_1_16800001_16860000", "DUP"),
    (GCNV_SEGMENTS, "CNV_1_10001_16800000", "CNV"),
    # gCNV joint cohort: DEL carrier plus no-call DUP carrier
    (GCNV_COHORT, "CNV_1_152555001_152590000", "CNV"),
    (GCNV_COHORT, "CNV_1_20100001_20150000", "DEL"),
    # DRAGEN legacy ASCN copy-neutral LOH, GT 1/2
    (DRAGEN_ASCN, "DRAGEN:CNLOH:1:152555002-152590000", "CNV"),
    # Canvas per-haplotype copy-number alleles
    (CANVAS, "Canvas:COMPLEXCNV:1:46001-111250", "CNV"),
    (CANVAS, "Canvas:LOSS:1:20100001-20150000", "CNV"),
    # single-ALT DRAGEN records keep the ALT type over SVTYPE=CNV
    (DRAGEN, "DRAGEN:LOSS:1:20100002-20150000", "DEL"),
    (DRAGEN, "DRAGEN:GAIN:1:16800002-16860000", "DUP"),
])
def test_record_type(vcf, record_id, expected):
    assert parse(vcf, record_id).event_type == expected


@pytest.mark.parametrize("vcf, record_id, gts", [
    (GCNV_SEGMENTS, "CNV_1_20100001_20150000", {"1": "DEL", "2": "DUP", "0": None, ".": None}),
    (GCNV_COHORT, "CNV_1_152555001_152590000",
     {"0/1": "DEL", "./2": "DUP", "1/2": "CNV", "./.": None, "0/0": None, "0|1": "DEL"}),
    (CANVAS, "Canvas:COMPLEXCNV:1:46001-111250", {"0/1": "CNV", "0/2": "CNV", "0/0": None}),
])
def test_type_for_gt(vcf, record_id, gts):
    variant = parse(vcf, record_id)
    assert {gt: variant.type_for_gt(gt) for gt in gts} == gts


@pytest.mark.parametrize("alt, expected", [
    ("<DEL>", "DEL"),
    ("<DUP:TANDEM>", "DUP"),
    ("<INS:ME:ALU>", "INS:ME:ALU"),
    ("<CNV>", "CNV"),
])
def test_single_alt_unchanged(alt, expected):
    line = f"1\t1000\tid\tN\t{alt}\t.\tPASS\tEND=5000\tGT\t0/1"
    variant = readVCFLine(line)
    assert variant.event_type == expected
    assert variant.type_for_gt("0/1") == expected


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    prefix = tmp_path_factory.mktemp("multi_alt") / "cnv"
    run("--build", "--files", GCNV_SEGMENTS, GCNV_COHORT, DRAGEN, DRAGEN_ASCN, CANVAS, GNOMAD,
        "--prefix", prefix)
    return Path(f"{prefix}.db")


def test_build_types_each_sample_by_gt(db):
    with sqlite3.connect(db) as con:
        rows = set(con.execute("SELECT sample, var, posA FROM SVDB"))
    assert rows == {
        # gCNV <= 4.1.4: reference segments (GT 0) are not stored
        ("HG002", "DUP", 16800001), ("HG002", "DEL", 20100001),
        # gCNV cohort: no-call DUP carrier (./.) is skipped, as for any ./.
        ("HG003", "DEL", 20100001), ("HG003", "DEL", 152555001),
        ("HG007", "DUP", 16800001), ("HG007", "DEL", 20100001),
        ("COLO829T", "CNV", 152555001),
        ("HG006", "CNV", 46000), ("NA12878", "CNV", 46000), ("HG006", "CNV", 20100000),
        ("gnomad_v2_1_mcnv", "CNV", 46000), ("gnomad_v2_1_mcnv", "CNV", 94000),
    }


def test_export_writes_cnv(db, tmp_path):
    run("--export", "--db", db, "--prefix", tmp_path / "out")
    text = (tmp_path / "out.vcf").read_text()
    assert "##ALT=<ID=CNV," in text
    assert ">,<" not in text
    cnv = [r for r in records(text) if r["INFO"]["SVTYPE"] == "CNV"]
    assert cnv and all(r["ALT"] == "<CNV>" for r in cnv)


def test_merge_gcnv_del_with_dragen_del():
    recs = records(run("--merge", "--vcf", DRAGEN, GCNV_SEGMENTS))
    at_del = [r for r in recs if r["POS"] == 20100001]
    assert len(at_del) == 1
    assert at_del[0]["INFO"]["FOUNDBY"] == "2"


def test_merge_mcnv_with_different_allele_counts():
    recs = records(run("--merge", "--vcf", GNOMAD, CANVAS))
    at_mcnv = [r for r in recs if r["POS"] == 46000]
    assert len(at_mcnv) == 1
    assert at_mcnv[0]["INFO"]["FOUNDBY"] == "2"


def test_query_db_dragen_del_hits_gcnv_del():
    recs = records(run("--query", "--db", GCNV_SEGMENTS, "--query_vcf", DRAGEN))
    occ = {r["POS"]: int(r["INFO"].get("OCC", 0)) for r in recs}
    assert occ[20100001] == 1


def test_query_sqdb_mcnv_hits_canvas(tmp_path):
    run("--build", "--files", GNOMAD, "--prefix", tmp_path / "gnomad")
    recs = records(run("--query", "--sqdb", tmp_path / "gnomad.db", "--query_vcf", CANVAS))
    occ = {r["POS"]: int(r["INFO"].get("OCC", 0)) for r in recs}
    assert occ[46000] == 1
