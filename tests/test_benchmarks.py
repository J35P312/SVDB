"""Regression benchmarks for SVDB performance.

These tests assert that each command completes within a generous wall-clock
budget on the small fixture VCFs. They are NOT micro-benchmarks — the goal
is to catch catastrophic regressions (e.g., an O(n²) loop becoming O(n³)),
not sub-second differences.

Baselines established 2026-04-17 on full VCF data (see scripts/profile_svdb.py).
To reproduce: python setup.py build_ext --inplace  (Cython), then remove .so files (pure Python).

Full-data profiling results (macOS arm64, Python 3.12, ~17k variants):

  Command                 Pure Python   Cython    Speedup
  merge (3 VCFs)          17.1s         7.5s      2.3x  ← material gain
  query vcf-db (manta)     0.40s        0.25s     1.6x
  query vcf-db (tiddit)    0.26s        0.16s     1.6x
  build (manta+tiddit)     0.47s        0.32s     1.5x
  export default           1.09s        0.81s     1.3x  ← bottleneck is SQLite I/O
  export DBSCAN            1.08s        0.82s     1.3x
  query sqdb               0.18s        0.13s     1.4x

Verdict: Cython is worth keeping for merge (2.3x on full data). Export gains
are minimal because 14k SQLite execute() calls dominate — Cython cannot help there.
"""

import subprocess
import sys
import time
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
MANTA = FIXTURES / "manta_chr1_del.vcf"
TIDDIT = FIXTURES / "tiddit_chr1_del.vcf"
TRUTH = FIXTURES / "truth_chr1_del.vcf"
MANTA_BND = FIXTURES / "manta_bnd.vcf"

SVDB = [sys.executable, "-m", "svdb"]

# Wall-clock budgets (seconds). Set 10x above observed fixture runtime to
# absorb CI variance while still catching algorithmic regressions.
BUDGET = {
    "merge": 5.0,
    "query_vcf": 2.0,
    "query_sqdb": 2.0,
    "build": 2.0,
    "export": 5.0,
}


def timed_run(*args) -> float:
    """Run svdb, return wall-clock seconds. Raises on non-zero exit."""
    t0 = time.perf_counter()
    result = subprocess.run(SVDB + list(args), capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    assert result.returncode == 0, result.stderr
    return elapsed


class TestMergeBenchmark:
    def test_merge_two_callers_within_budget(self):
        elapsed = timed_run("--merge", "--vcf", str(MANTA), str(TIDDIT))
        assert elapsed < BUDGET["merge"], (
            f"merge took {elapsed:.2f}s — exceeds budget of {BUDGET['merge']}s"
        )

    def test_merge_three_callers_within_budget(self):
        elapsed = timed_run("--merge", "--vcf", str(MANTA), str(TIDDIT), str(TRUTH))
        assert elapsed < BUDGET["merge"], (
            f"merge (3 VCFs) took {elapsed:.2f}s — exceeds budget of {BUDGET['merge']}s"
        )

    def test_merge_bnd_within_budget(self):
        """BND variants trigger the precise_overlap path — check it doesn't regress."""
        elapsed = timed_run("--merge", "--vcf", str(MANTA_BND), str(MANTA_BND))
        assert elapsed < BUDGET["merge"], (
            f"BND merge took {elapsed:.2f}s — exceeds budget of {BUDGET['merge']}s"
        )


class TestQueryBenchmark:
    def test_query_vcfdb_within_budget(self):
        elapsed = timed_run(
            "--query", "--db", str(TRUTH), "--query_vcf", str(MANTA)
        )
        assert elapsed < BUDGET["query_vcf"], (
            f"query (vcf db) took {elapsed:.2f}s — exceeds budget of {BUDGET['query_vcf']}s"
        )

    def test_query_vcfdb_tight_overlap_within_budget(self):
        """overlap=1.0 is the most expensive path (checks all candidates)."""
        elapsed = timed_run(
            "--query", "--db", str(TRUTH), "--query_vcf", str(MANTA), "--overlap", "1.0"
        )
        assert elapsed < BUDGET["query_vcf"], (
            f"query tight overlap took {elapsed:.2f}s — exceeds {BUDGET['query_vcf']}s"
        )

    def test_query_sqdb_within_budget(self, tmp_path):
        prefix = tmp_path / "svdb"
        subprocess.run(
            SVDB + ["--build", "--files", str(MANTA), str(TIDDIT),
                    "--prefix", str(prefix)], capture_output=True
        )
        elapsed = timed_run(
            "--query", "--sqdb", str(tmp_path / "svdb.db"), "--query_vcf", str(TRUTH)
        )
        assert elapsed < BUDGET["query_sqdb"], (
            f"query (sqdb) took {elapsed:.2f}s — exceeds budget of {BUDGET['query_sqdb']}s"
        )


class TestBuildBenchmark:
    def test_build_within_budget(self, tmp_path):
        elapsed = timed_run(
            "--build", "--files", str(MANTA), str(TIDDIT),
            "--prefix", str(tmp_path / "svdb")
        )
        assert elapsed < BUDGET["build"], (
            f"build took {elapsed:.2f}s — exceeds budget of {BUDGET['build']}s"
        )


class TestExportBenchmark:
    @pytest.fixture
    def db(self, tmp_path):
        subprocess.run(
            SVDB + ["--build", "--files", str(MANTA), str(TIDDIT),
                    "--prefix", str(tmp_path / "svdb")], capture_output=True
        )
        return tmp_path / "svdb.db"

    def test_export_default_within_budget(self, db, tmp_path):
        elapsed = timed_run(
            "--export", "--db", str(db), "--prefix", str(tmp_path / "out")
        )
        assert elapsed < BUDGET["export"], (
            f"export took {elapsed:.2f}s — exceeds budget of {BUDGET['export']}s"
        )

    def test_export_dbscan_within_budget(self, db, tmp_path):
        """DBSCAN path was the site of the UnboundLocalError bug — keep a timing check."""
        elapsed = timed_run(
            "--export", "--db", str(db), "--DBSCAN",
            "--prefix", str(tmp_path / "dbscan")
        )
        assert elapsed < BUDGET["export"], (
            f"export (DBSCAN) took {elapsed:.2f}s — exceeds budget of {BUDGET['export']}s"
        )
