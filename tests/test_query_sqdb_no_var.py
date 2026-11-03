"""Integration tests for issue #108: --no_var must apply to --query --sqdb as it does to --db.

Fixtures (tests/fixtures/sqdb_no_var): one sample with a DUP on chr1 and a DEL on chr2;
the query has a DEL at the same position on each chromosome.
"""

import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "sqdb_no_var"
SAMPLE = FIXTURES / "sample.vcf"
QUERY = FIXTURES / "query.vcf"

SVDB = [sys.executable, "-m", "svdb"]


def run(*args):
    r = subprocess.run(SVDB + list(args), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def occ_by_id(vcf_text):
    """Map record ID -> OCC value (None if not annotated)."""
    out = {}
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        out[f[2]] = info.get("OCC")
    return out


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    prefix = tmp_path_factory.mktemp("sqdb") / "svdb"
    run("--build", "--files", str(SAMPLE), "--prefix", str(prefix))
    return prefix.with_suffix(".db")


def query_sqdb(db, *opts):
    return occ_by_id(run("--query", "--query_vcf", str(QUERY), "--sqdb", str(db), *opts))


def test_default_matches_same_type_only(db):
    assert query_sqdb(db) == {"del_vs_dup": None, "del_vs_del": "1"}


def test_no_var_matches_across_types(db):
    assert query_sqdb(db, "--no_var") == {"del_vs_dup": "1", "del_vs_del": "1"}


def test_no_var_sqdb_matches_vcf_db(db, tmp_path):
    exported = tmp_path / "exported"
    run("--export", "--db", str(db), "--prefix", str(exported))
    vcf_db = occ_by_id(run("--query", "--query_vcf", str(QUERY), "--db", f"{exported}.vcf",
                           "--in_occ", "OCC", "--in_frq", "FRQ", "--no_var"))
    assert query_sqdb(db, "--no_var") == vcf_db
