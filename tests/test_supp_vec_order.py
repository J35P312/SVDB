"""SUPP_VEC digits follow the --vcf order, as labelled in the header (#116).

Fixtures (tests/fixtures/supp_vec_order): zeta.vcf and alpha.vcf, so that
--vcf order differs from alphabetical tag order. Each file has one INS of its
own (POS 1000000 in zeta, 2000000 in alpha) and one INS shared at 3000000.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "supp_vec_order"
ZETA = FIXTURES / "zeta.vcf"
ALPHA = FIXTURES / "alpha.vcf"

# POS -> tags whose file holds the variant
EXPECTED = {"1000000": {"zeta"}, "2000000": {"alpha"}, "3000000": {"zeta", "alpha"}}


def merge(*args):
    r = subprocess.run([sys.executable, "-m", "svdb", "--merge", *map(str, args)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def header_order(vcf_text):
    m = re.search(r"ID=SUPP_VEC,.*\(order: ([^)]*)\)", vcf_text)
    assert m, "no SUPP_VEC header line"
    return m.group(1).split()


def supp_vecs(vcf_text):
    """POS -> SUPP_VEC string."""
    out = {}
    for line in vcf_text.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        out[f[1]] = re.search(r"(?:^|;)SUPP_VEC=([01]+)", f[7]).group(1)
    return out


CASES = {
    "priority_same_as_vcf": ([f"{ZETA}:zeta", f"{ALPHA}:alpha"], ["--priority", "zeta,alpha"], ["zeta", "alpha"]),
    "priority_reversed": ([f"{ZETA}:zeta", f"{ALPHA}:alpha"], ["--priority", "alpha,zeta"], ["zeta", "alpha"]),
    "no_priority": ([ZETA, ALPHA], [], ["zeta", "alpha"]),
    "vcf_alphabetical": ([f"{ALPHA}:alpha", f"{ZETA}:zeta"], ["--priority", "zeta,alpha"], ["alpha", "zeta"]),
}


@pytest.mark.parametrize("vcfs,extra,vcf_order", CASES.values(), ids=CASES.keys())
def test_supp_vec_matches_header(vcfs, extra, vcf_order):
    out = merge("--vcf", *vcfs, *extra)
    order = header_order(out)
    assert order == vcf_order
    vecs = supp_vecs(out)
    assert set(vecs) == set(EXPECTED)
    for pos, tags in EXPECTED.items():
        assert vecs[pos] == "".join("1" if t in tags else "0" for t in order), pos


def test_supp_vec_independent_of_priority():
    vcfs = [f"{ZETA}:zeta", f"{ALPHA}:alpha"]
    a = merge("--vcf", *vcfs, "--priority", "zeta,alpha")
    b = merge("--vcf", *vcfs, "--priority", "alpha,zeta")
    assert header_order(a) == header_order(b)
    assert supp_vecs(a) == supp_vecs(b)
