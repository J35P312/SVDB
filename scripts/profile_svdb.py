"""Profile SVDB commands on real VCF data to establish a pure-Python baseline.

Run with:
    python scripts/profile_svdb.py [vcf_dir]

vcf_dir defaults to ../vcf_files relative to repo root.
Output: profiling stats per command, sorted by cumulative time.

Purpose: compare against a Cython-compiled run once setup.py compilation
is confirmed working. The Cython modules are:
  build_module, overlap_module, dbscan, read_vcf,
  merge_vcf_module_cython, query_module, export_module
"""

import argparse
import cProfile
import io
import pstats
import sys
import tempfile
from pathlib import Path


def _run_profiled(argv: list[str]) -> tuple[pstats.Stats, io.StringIO]:
    import importlib
    import svdb.__main__ as entry
    importlib.reload(entry)

    sys.argv = argv
    pr = cProfile.Profile()
    pr.enable()
    # suppress VCF stdout during profiling
    _real_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        entry.main()
    except SystemExit:
        pass
    finally:
        sys.stdout = _real_stdout
    pr.disable()

    stream = io.StringIO()
    stats = pstats.Stats(pr, stream=stream).sort_stats("cumulative")
    return stats, stream


def print_profile(label: str, stats: pstats.Stats, stream: io.StringIO, n: int = 20):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    stats.print_stats(n)
    print(stream.getvalue())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf_dir", nargs="?", default="../vcf_files",
                        help="Directory containing manta.vcf, NA12878.tiddit.pass.vcf, "
                             "Personalis_1000_Genomes_deduplicated_deletions.vcf")
    parser.add_argument("--top", type=int, default=15,
                        help="Number of top functions to show per command (default: 15)")
    args = parser.parse_args()

    vcf_dir = Path(args.vcf_dir)
    manta = vcf_dir / "manta.vcf"
    tiddit = vcf_dir / "NA12878.tiddit.pass.vcf"
    truth = vcf_dir / "Personalis_1000_Genomes_deduplicated_deletions.vcf"

    for f in [manta, tiddit, truth]:
        if not f.exists():
            print(f"ERROR: {f} not found — pass the correct vcf_dir", file=sys.stderr)
            sys.exit(1)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        commands = [
            ("merge (3 VCFs)",
             ["svdb", "--merge", "--vcf", str(manta), str(tiddit), str(truth)]),

            ("query vcf-db (tiddit vs truth)",
             ["svdb", "--query", "--db", str(truth), "--query_vcf", str(tiddit)]),

            ("query vcf-db (manta vs truth)",
             ["svdb", "--query", "--db", str(truth), "--query_vcf", str(manta)]),

            ("build (manta + tiddit)",
             ["svdb", "--build", "--files", str(manta), str(tiddit),
              "--prefix", str(tmp / "svdb")]),

            ("export default",
             ["svdb", "--export", "--db", str(tmp / "svdb.db"),
              "--prefix", str(tmp / "export")]),

            ("export DBSCAN",
             ["svdb", "--export", "--db", str(tmp / "svdb.db"),
              "--DBSCAN", "--prefix", str(tmp / "export_dbscan")]),

            ("query sqdb (truth vs manta+tiddit db)",
             ["svdb", "--query", "--sqdb", str(tmp / "svdb.db"),
              "--query_vcf", str(truth)]),
        ]

        # build must run before export/sqdb query
        results = {}
        for label, argv in commands:
            print(f"Profiling: {label} ...", file=sys.stderr)
            stats, stream = _run_profiled(argv)
            results[label] = (stats, stream)

        for label, (stats, stream) in results.items():
            print_profile(label, stats, stream, n=args.top)


if __name__ == "__main__":
    main()
