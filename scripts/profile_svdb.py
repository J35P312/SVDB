"""Profile SVDB commands on real VCF data.

Usage
-----
Copy scripts/profile_config.toml.example to scripts/profile_config.toml,
fill in your file paths and caller names, then run:

    python scripts/profile_svdb.py [--config PATH] [--top N] [--sort KEY]

Options
-------
--config PATH   TOML config file describing the VCFs to profile.
                Defaults to scripts/profile_config.toml next to this script.
--top N         Number of top functions to show per command (default: 15).
--sort KEY      cProfile sort key: cumulative (default), tottime, calls.

Config format
-------------
See profile_config.toml.example for the full schema. In brief:

  [vcfs]
  caller_name = "/abs/path/to/caller.vcf"
  ...

  [options]
  truth = "caller_name"          # used as --db for vcf-db query commands
  merge_keys = ["a", "b", ...]   # subset to merge (default: all non-truth)
  build_keys  = ["a", "b", ...]  # subset to build  (default: all non-truth)
  no_var_merge = false            # pass --no_var to the multi-VCF merge

Purpose
-------
Catches algorithmic regressions and validates optimisation work.  Always
profiles the local checkout (repo root is inserted into sys.path before
any svdb import).  Wall-clock numbers include cProfile overhead; ratios
between runs on the same machine are reliable, absolute numbers are not.
"""

import argparse
import cProfile
import io
import pstats
import sys
import tempfile
from pathlib import Path

# Always profile the local checkout, not any installed package.
_REPO_ROOT = str(Path(__file__).parent.parent.resolve())
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

try:
    import tomllib  # stdlib on Python 3.11+
except ImportError:
    try:
        import tomli as tomllib  # pip install tomli  (Python 3.9 / 3.10)
    except ImportError:
        tomllib = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Profiling helpers
# ---------------------------------------------------------------------------

def _run_profiled(argv: list[str]) -> tuple[pstats.Stats, io.StringIO]:
    """Run an svdb command in-process under cProfile, return stats."""
    import importlib
    import svdb.__main__ as entry
    importlib.reload(entry)

    sys.argv = argv
    pr = cProfile.Profile()
    pr.enable()
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


def _print_profile(
    label: str, stats: pstats.Stats, stream: io.StringIO,
    n: int, sort_key: str,
) -> None:
    print(f"\n{'=' * 60}")
    print(f"  {label}")
    print(f"{'=' * 60}")
    stats.sort_stats(sort_key)
    stats.print_stats(n)
    print(stream.getvalue())


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def _load_config(path: Path) -> dict:
    if tomllib is None:
        print(
            "ERROR: TOML support not available.\n"
            "  Python 3.11+ has tomllib built-in.\n"
            "  For Python 3.9/3.10: pip install tomli",
            file=sys.stderr,
        )
        sys.exit(1)
    if not path.exists():
        print(
            f"ERROR: config file not found: {path}\n"
            f"  Copy scripts/profile_config.toml.example to {path} and fill in your paths.",
            file=sys.stderr,
        )
        sys.exit(1)
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def _resolve_vcfs(config: dict) -> dict[str, Path]:
    """Return {name: Path} for every VCF in [vcfs], erroring on missing files."""
    raw = config.get("vcfs", {})
    if not raw:
        print("ERROR: [vcfs] section is empty in config.", file=sys.stderr)
        sys.exit(1)
    vcfs: dict[str, Path] = {}
    missing = []
    for name, path_str in raw.items():
        p = Path(path_str)
        if not p.exists():
            missing.append(f"  {name}: {path_str}")
        else:
            vcfs[name] = p
    if missing:
        print("ERROR: the following VCF files were not found:", file=sys.stderr)
        print("\n".join(missing), file=sys.stderr)
        sys.exit(1)
    return vcfs


# ---------------------------------------------------------------------------
# Command generation
# ---------------------------------------------------------------------------

def _build_commands(vcfs: dict[str, Path], opts: dict, tmp: Path) -> list[tuple[str, list[str]]]:
    truth_key = opts.get("truth")
    merge_keys = opts.get("merge_keys") or [k for k in vcfs if k != truth_key]
    build_keys = opts.get("build_keys") or [k for k in vcfs if k != truth_key]
    no_var = opts.get("no_var_merge", False)

    merge_vcfs = [str(vcfs[k]) for k in merge_keys if k in vcfs]
    build_vcfs = [str(vcfs[k]) for k in build_keys if k in vcfs]
    truth_vcf = str(vcfs[truth_key]) if truth_key and truth_key in vcfs else None

    commands: list[tuple[str, list[str]]] = []

    # --- merge ---
    if len(merge_vcfs) >= 2:
        label = f"merge ({len(merge_vcfs)} VCFs: {', '.join(merge_keys)})"
        cmd = ["svdb", "--merge", "--vcf"] + merge_vcfs
        if no_var:
            cmd.insert(2, "--no_var")
        commands.append((label, cmd))

    # --- single large file merged with itself (stress-test) ---
    # Pick the build VCF with the most variants if we can stat file size as proxy
    if build_vcfs:
        largest = max(build_vcfs, key=lambda p: Path(p).stat().st_size)
        largest_name = next(k for k in build_keys if str(vcfs[k]) == largest)
        commands.append((
            f"merge ({largest_name} × 2, self-merge stress)",
            ["svdb", "--merge", "--vcf", largest, largest],
        ))

    # --- build ---
    if build_vcfs:
        label = f"build ({len(build_vcfs)} VCFs: {', '.join(build_keys)})"
        commands.append((label, [
            "svdb", "--build", "--files"] + build_vcfs + ["--prefix", str(tmp / "svdb")],
        ))
        commands.append(("export (default)", [
            "svdb", "--export", "--db", str(tmp / "svdb.db"),
            "--prefix", str(tmp / "export"),
        ]))
        commands.append(("export (DBSCAN)", [
            "svdb", "--export", "--db", str(tmp / "svdb.db"),
            "--DBSCAN", "--prefix", str(tmp / "export_dbscan"),
        ]))

    # --- query ---
    if truth_vcf and build_vcfs:
        commands.append((
            f"query vcf-db (truth={truth_key} vs {build_keys[0]})",
            ["svdb", "--query", "--db", str(vcfs[build_keys[0]]),
             "--query_vcf", truth_vcf],
        ))
        if (tmp / "svdb.db").exists():
            commands.append((
                f"query sqdb (truth={truth_key} vs built db)",
                ["svdb", "--query", "--sqdb", str(tmp / "svdb.db"),
                 "--query_vcf", truth_vcf],
            ))

    return commands


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    _default_config = Path(__file__).parent / "profile_config.toml"

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config", type=Path, default=_default_config,
        help=f"TOML config file (default: {_default_config})",
    )
    parser.add_argument(
        "--top", type=int, default=15,
        help="Number of top functions to show per command (default: 15)",
    )
    parser.add_argument(
        "--sort", default="cumulative",
        choices=["cumulative", "tottime", "calls"],
        help="cProfile sort key (default: cumulative)",
    )
    args = parser.parse_args()

    config = _load_config(args.config)
    vcfs = _resolve_vcfs(config)
    opts = config.get("options", {})

    print(f"Loaded {len(vcfs)} VCF(s) from {args.config}:", file=sys.stderr)
    for name, path in vcfs.items():
        print(f"  {name}: {path}", file=sys.stderr)

    with tempfile.TemporaryDirectory() as _tmp:
        tmp = Path(_tmp)
        commands = _build_commands(vcfs, opts, tmp)

        # build must run before export/sqdb-query; execute in order, collect results
        results: dict[str, tuple[pstats.Stats, io.StringIO]] = {}
        for label, argv in commands:
            print(f"Profiling: {label} ...", file=sys.stderr)
            stats, stream = _run_profiled(argv)
            results[label] = (stats, stream)

        for label, (stats, stream) in results.items():
            _print_profile(label, stats, stream, n=args.top, sort_key=args.sort)


if __name__ == "__main__":
    main()
