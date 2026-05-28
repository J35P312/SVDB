import argparse
import importlib.metadata
import logging
import os
import sys

from . import build_module, export_module, merge_vcf_module, query_module

logger = logging.getLogger(__name__)


def _fraction(flag: str):
    def _parse(value: str) -> float:
        f = float(value)
        if not 0.0 <= f <= 1.0:
            raise argparse.ArgumentTypeError(f"{flag} must be in [0.0, 1.0], got {f}")
        return f
    return _parse


def _positive_int(flag: str):
    def _parse(value: str) -> int:
        try:
            n = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{flag} must be a whole number ≥ 1, got {value!r}")
        if n < 1:
            raise argparse.ArgumentTypeError(f"{flag} must be ≥ 1, got {n}")
        return n
    return _parse


def _max_ins_seq_len(value: str):
    """Parse --max_ins_seq_len: positive integer, or 0 meaning unlimited (None)."""
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--max_ins_seq_len must be a positive integer or 0 (unlimited), got {value!r}"
        )
    if n < 0:
        raise argparse.ArgumentTypeError(
            f"--max_ins_seq_len must be ≥ 0 (0 = unlimited), got {n}"
        )
    return None if n == 0 else n


def _setup_logging(debug: bool) -> None:
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )


def _add_ins_flags(parser: argparse.ArgumentParser) -> None:
    """Add the shared insertion matching flags to a subcommand parser."""
    parser.add_argument(
        '--ins_distance', type=int, default=None,
        help="maximum distance to match two insertions "
             "(default: 25; profile cohort/position_only: 50)")
    parser.add_argument(
        '--ins_svlen_ratio', type=_fraction('--ins_svlen_ratio'), default=None,
        help="minimum SVLEN ratio (min/max) for insertions with known length "
             "(default: 0.90; profile cohort: 0.80)")
    parser.add_argument(
        '--ins_seq_similarity', type=_fraction('--ins_seq_similarity'), default=None,
        help="minimum Levenshtein sequence similarity (0–1); explicit value "
             "overrides --data_profile for this parameter (effective default: 0.75)")
    parser.add_argument(
        '--data_profile', choices=["sample", "cohort", "position_only"], default=None,
        help="insertion matching preset. "
             "sample: strict (dist=25, ratio=0.90, sim=0.85) — same individual / technology. "
             "cohort: permissive (dist=50, ratio=0.80, sim=0.75) — cross-individual or cross-caller. "
             "position_only: no sequence gate, match on position+SVLEN only (dist=50, ratio=0.90). "
             "Individual --ins_* flags override profile values.")
    parser.add_argument(
        '--max_ins_seq_len', type=_max_ins_seq_len, default=1000,
        help="sequences longer than N bp are excluded from sequence similarity "
             "and fall back to position+SVLEN (default: 1000); "
             "use 0 for no cap (compare all sequences regardless of length).")


def make_query_calls (args, queries, keyword):
    if len(queries) > 1 and args.prefix:
        if all(variable is not None for variable in [args.in_occ, args.out_occ, args.in_frq, args.out_frq]):
            in_occs     = args.in_occ.split(",")
            in_frqs     = args.in_frq.split(",")
            out_occs    = args.out_occ.split(",")
            out_frqs    = args.out_frq.split(",")
            orig_prefix = args.prefix
            if (len(queries) == len(in_occs) == len(in_frqs) == len(out_occs) == len(out_frqs)):
                for ind in range(len(queries)):
                    if keyword   == "db":
                        args.db      = queries[ind]
                    elif keyword == "sqdb":
                        args.sqdb    = queries[ind]
                    elif keyword == "bedpedb":
                        args.bedpedb = queries[ind]
                    args.in_occ      = None if in_occs[ind] == "default" else in_occs[ind]
                    args.in_frq      = None if in_frqs[ind] == "default" else in_frqs[ind]
                    args.out_occ     = out_occs[ind]
                    args.out_frq     = out_frqs[ind]
                    if ind < len(queries)-1:
                        args.prefix  = orig_prefix + "_" + str(ind)
                        output_file  = args.prefix + "_query.vcf_tmp"
                    else:
                        args.prefix = orig_prefix
                        output_file  = args.prefix + "_query.vcf"
                    query_module.main(args, output_file)
                    if ind > 0:
                        os.remove(args.query_vcf)
                    args.query_vcf = output_file
        else:
            logger.error("please ensure that both count and frequency tags are specified for all samples")
            sys.exit(1)
    elif len(queries) > 1 and not args.prefix:
        logger.error("please provide a --prefix when querying multiple databases")
        sys.exit(1)
    elif len(queries) == 1 and args.prefix:
        output_file  = args.prefix + "_query.vcf"
        query_module.main(args, output_file)
    else:
        query_module.main(args)

def main():
    version = importlib.metadata.version("svdb")
    parser = argparse.ArgumentParser(
        f"""SVDB-{version}: use --build to construct databases, --query to annotate a VCF, --merge to merge callers, or --export to export a database to VCF""", add_help=False)
    parser.add_argument('--build', help="create a db",
                        required=False, action="store_true")
    parser.add_argument('--query', help="query a db",
                        required=False, action="store_true")
    parser.add_argument('--merge', help="merge similar structural variants within a vcf file",
                        required=False, action="store_true")
    parser.add_argument('--export', help="export a database",
                        required=False, action="store_true")
    parser.add_argument('--debug', help="enable debug logging to stderr",
                        required=False, action="store_true")
    args, unknown = parser.parse_known_args()
    _setup_logging(args.debug)

    if args.query:
        parser = argparse.ArgumentParser(
            f"""SVDB-{version}: query module""")
        parser.add_argument('--query', help="query a db", required=False, action="store_true")
        parser.add_argument('--query_vcf', type=str, help="a vcf used to query the db", required=True)
        parser.add_argument('--db', type=str, help="path to a db vcf, or a comma-separated list of vcfs")
        parser.add_argument('--sqdb', type=str, help="path to a SVDB sqlite db, or a comma-separated list of dbs")
        parser.add_argument('--bedpedb', type=str,
                            help="path to a SV db in chrA-posA-chrB-posB-type-count-frequency format, or a comma-separated list of files")
        parser.add_argument('--in_occ', type=str,
                            help="the allele count tag in the db INFO column (usually OCC or AN); "
                                 "required when querying multiple databases. Use 'default' to use the default tag for a specific db")
        parser.add_argument('--in_frq', type=str,
                            help="the allele frequency tag in the db INFO column (usually FRQ or AF); "
                                 "required when querying multiple databases. Use 'default' to use the default tag for a specific db")
        parser.add_argument('--out_occ', type=str, default="OCC",
                            help="output tag for allele count (default: OCC); required when querying multiple databases")
        parser.add_argument('--out_frq', type=str, default="FRQ",
                            help="output tag for allele frequency (default: FRQ); required when querying multiple databases")
        parser.add_argument('--max_frq', type=_fraction('--max_frq'), default=1,
                            help='only include variants with frequency at or below this value (default: 1, i.e. all variants)')
        parser.add_argument('--prefix', type=str, default=None,
                            help="prefix for the output file (default: print to stdout); required when querying multiple databases")
        parser.add_argument('--bnd_distance', type=int, default=10000,
                            help="maximum distance between two similar breakpoints (default: 10000)")
        _add_ins_flags(parser)
        parser.add_argument('--overlap', type=_fraction('--overlap'), default=0.6,
                            help="minimum reciprocal overlap required to match two events "
                                 "(0 = anything touching; 1 = identical) (default: 0.6)")
        parser.add_argument('--memory',
                            help="load the database into memory: increases memory use but lowers query time (sqdb only)", required=False, action="store_true")
        parser.add_argument('--no_var',
                            help="count overlapping variants of different type as hits in the db", required=False, action="store_true")
        parser.add_argument('--debug', help=argparse.SUPPRESS,
                            required=False, action="store_true")
        args = parser.parse_args()
        args.version = version

        if(args.db or args.sqdb or args.bedpedb):
            if(args.db):
                queries = args.db.split(",")
                make_query_calls(args, queries, "db")
            if(args.sqdb):
                queries = args.sqdb.split(",")
                make_query_calls(args, queries, "sqdb")
            if(args.bedpedb):
                queries = args.bedpedb.split(",")
                make_query_calls(args, queries, "bedpedb")
        else:
            logger.error("invalid db option, choose --db to use the vcf db or sqdb to use the sqlite db")

    elif args.build:
        parser = argparse.ArgumentParser(
            f"""SVDB-{version}: build module""")
        parser.add_argument('--build', help="create a db",
                            required=False, action="store_true")
        parser.add_argument(
            '--pass_only', help="only include variants with PASS or . in the FILTER field", required=False, action="store_true")
        parser.add_argument(
            '--passonly', help=argparse.SUPPRESS, required=False, action="store_true")
        parser.add_argument('--files', type=str, nargs='*',
                            help="create a db using the specified vcf files (cannot be used with --folder)")
        parser.add_argument(
            '--folder', type=str, help="create a db using all vcf files in the given folder")
        parser.add_argument('--prefix', type=str, default="SVDB",
                            help="prefix for the output file (default: SVDB)")
        parser.add_argument('--upgrade', help="create the INS table and backfill insertion sequences "
                            "from the provided VCFs; schema-only, no existing SVDB rows are changed. "
                            "Exits with INFO if the INS table already exists. "
                            "Warns for DB samples with no matching VCF; logs DB-absent VCF samples.",
                            required=False, action="store_true")
        parser.add_argument('--max_ins_seq_len', type=int, default=None,
                            help="maximum insertion sequence length (bp) to store; sequences longer than this "
                                 "are stored with NULL sequence but retain SVLEN for length-ratio matching "
                                 "(default: no limit)")
        parser.add_argument('--debug', help=argparse.SUPPRESS,
                            required=False, action="store_true")
        args = parser.parse_args()
        if args.passonly:
            logger.warning("--passonly is deprecated; use --pass_only instead")
            args.pass_only = True
        args.version = version
        if (args.files and args.folder):
            logger.error("only one DB build input source may be selected (--files or --folder)")
            sys.exit(1)

        if args.upgrade and not args.files and not args.folder:
            logger.error("--upgrade requires --files or --folder to backfill insertion sequences")
            sys.exit(1)
        if args.upgrade:
            build_module.main(args)
        elif args.folder or args.files:
            build_module.main(args)
        else:
            logger.error("use --files or --folder to provide input for the database creation algorithm")
    elif args.export:
        parser = argparse.ArgumentParser(
            f"""SVDB-{version}: export module; export the variants of the SVDB sqlite database into a vcf file""")
        parser.add_argument('--export', help="export the database",
                            required=False, action="store_true")
        parser.add_argument('--db', type=str, required=True,
                            help="the SQLite database")
        parser.add_argument(
            '--no_merge', help="skip the merging of variants, print all variants in the db to a vcf file", required=False, action="store_true")
        parser.add_argument('--bnd_distance', type=int, default=2500,
                            help="maximum distance between two similar precise breakpoints (default: 2500)")
        _add_ins_flags(parser)
        parser.add_argument('--overlap', type=_fraction('--overlap'), default=0.8,
                            help="minimum reciprocal overlap required to merge two events "
                                 "(0 = anything touching; 1 = identical) (default: 0.8)")
        parser.add_argument(
            '--coarse', help="skip the second-pass refinement (overlap/SVLEN/sequence gates) and use "
                             "centroid-based representative selection directly from the DBSCAN spatial "
                             "clusters; produces fewer, coarser clusters. Combine with --epsilon/--min_pts to tune.",
            required=False, action="store_true")
        parser.add_argument(
            '--DBSCAN', help=argparse.SUPPRESS, required=False, action="store_true")
        parser.add_argument('--epsilon', type=float, default=500,
                            help="used together with --coarse; spatial grouping radius in bp (DBSCAN-style "
                                 "epsilon): variants within this distance of each other are candidates "
                                 "for the same cluster (default: 500)", required=False)
        parser.add_argument('--min_pts', type=_positive_int('--min_pts'), default=2,
                            help="used together with --coarse; DBSCAN-style min_pts: minimum number of "
                                 "consecutively-positioned variants that must all fall within --epsilon "
                                 "to seed a cluster — isolated variants become singletons "
                                 "(default: 2, meaning any pair within --epsilon forms a cluster)", required=False)
        parser.add_argument('--prefix', type=str, default="SVDB",
                            help="prefix for the output file (default: same as input)")
        parser.add_argument(
            '--memory', help="load the database into memory: increases memory use but lowers export time", required=False, action="store_true")
        parser.add_argument('--strip_chr', help="strip the 'chr' prefix from chromosome names in the output VCF",
                            required=False, action="store_true")
        parser.add_argument('--samples', choices=['on', 'off'], default='on',
                            help="include per-sample genotype columns (default: on); use 'off' for sites-only output analogous to gnomAD --sites-only")
        parser.add_argument('--cluster_method', choices=['star', 'union_find'], default='star',
                            help="second-pass clustering algorithm applied within each DBSCAN spatial group: "
                                 "'star' = greedy star, highest-degree representative, no transitivity (default); "
                                 "'union_find' = transitive closure, fewer output clusters, higher OCC counts")
        parser.add_argument('--workers', type=int, default=0,
                            help="parallel worker processes for clustering (default: 0 = use all logical CPUs; 1 = serial). "
                                 "To find your optimal N: time with --workers 1, then increase until wall-clock time stops improving — "
                                 "that is the serial floor (DB fetch, DBSCAN, I/O). On shared systems set N explicitly to be a good neighbour.")
        parser.add_argument('--debug', help=argparse.SUPPRESS,
                            required=False, action="store_true")
        args = parser.parse_args()

        if args.DBSCAN:
            logger.warning("--DBSCAN is deprecated; use --coarse instead")
            args.coarse = True

        # merging will be impossible
        if args.no_merge:
            args.overlap = float("inf")
            args.bnd_distance = -1

        args.version = version
        export_module.main(args)

    elif args.merge:
        parser = argparse.ArgumentParser(
            f"""SVDB-{version}: vcf_merge module""")
        parser.add_argument(
            '--merge', help="merge structural variants", required=False, action="store_true")
        parser.add_argument(
            '--no_tag', help="do not add the VARID and set entries to the INFO field", required=False, action="store_true")
        parser.add_argument(
            '--notag', help=argparse.SUPPRESS, required=False, action="store_true")
        parser.add_argument('--vcf', nargs='*', type=str,
                            help="input vcf files to merge; use --priority to set caller precedence", required=True)
        parser.add_argument(
            '--priority', type=str, help="prioritise the input files, using the format --vcf caller1.vcf:2 caller2.vcf:1 --priority 1,2")
        parser.add_argument('--bnd_distance', type=int, default=2000,
                            help="maximum distance between two similar precise breakpoints (default: 2000)")
        _add_ins_flags(parser)
        parser.add_argument('--overlap', type=_fraction('--overlap'), default=0.95,
                            help="minimum reciprocal overlap required to merge two events "
                                 "(0 = anything touching; 1 = identical) (default: 0.95)")
        parser.add_argument(
            '--no_intra', help="no merging of variants within the same vcf", required=False, action="store_true")
        parser.add_argument(
            '--no_var', help="variants of different type will be merged", required=False, action="store_true")
        parser.add_argument(
            '--pass_only', help="merge only variants labeled PASS", required=False, action="store_true")
        parser.add_argument(
            '--same_order', help="across all input vcf files, the order of the sample columns is the same", required=False, action="store_true")
        parser.add_argument('--debug', help=argparse.SUPPRESS,
                            required=False, action="store_true")
        args = parser.parse_args()
        args.version = version
        if args.notag:
            logger.warning("--notag is deprecated; use --no_tag instead")
            args.no_tag = True
        result = merge_vcf_module.main(args)
        if result is not None and result < 0:
            sys.exit(1)

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
