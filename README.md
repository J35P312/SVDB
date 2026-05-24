# SVDB

SVDB is a toolkit for constructing and querying structural variant databases. The databases are constructed using the output vcf files from structural variant callers such as TIDDIT, Manta, Fermikit or Delly.
SVDB may also be used to merge SV vcf files from multiple callers or individuals.

## Supported public databases

SVDB query supports public databases such as thousand genomes SV map and Gnomad SV, as well as most multisample SV vcf files.

The thousand genomes SV database:
<ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/>

The swegen SVDB:
<https://swefreq.nbis.se/>

The GNOMAD SV database:
<https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2_sv.sites.vcf.gz>

External databases are run like this:

```bash
svdb --query \
     --query_vcf /home/jesper/vcf/6_pair_limit/P2109_120.clean.dedup.recal_FindSV.vcf \
     --out_occ GNOMAD_AC \
     --out_frq GNOMAD_AF \
     --in_occ AN \
     --in_frq AF \
     --db /home/jesper/Downloads/gnomad_sv/gnomad_v2_sv.sites.vcf
```

Here the `AF` and `AN` are the allele frequency tags of the database, the `AF` is a float, and `AN` is an integer. These tags will be added to the annotated output vcf, and named `GNOMAD_AC`, `GNOMAD_AF`.

## Install

Dependencies: SVDB requires Python 3.9+ and numpy.

Install from PyPI:

```bash
pip install svdb
```

Install from source:

```bash
git clone https://github.com/J35P312/SVDB.git
cd SVDB
pip install .
```

SVDB is available as a container on BioContainers:

<https://quay.io/repository/biocontainers/svdb?tab=tags&tag=>

SVDB is available on Singularity:

```bash
singularity pull shub://J35P312/SVDB
```

## Modules

SVDB consists of modules that are used to build, query, export, and analyse structural variant databases. These are the modules:

## Build

This module is used to construct structural variant databases from vcf files. The database may then be queried to compute the frequency of structural variants, or exported into a vcf file. These are the commands used to construct a structural variation database:

```text
print a help message
    svdb  --build --help
Construct a database from a set of vcf files:
    svdb --build --files sample1.vcf sample2.vcf sample3.vcf
Construct a database from vcf files stored in a folder:
    svdb --build --folder SV_analysis_folder/
Upgrade an existing database schema to the current SVDB version:
    svdb --build --upgrade --prefix existing_db
Upgrade schema and backfill insertion data from the original VCFs:
    svdb --build --upgrade --files sample1.vcf sample2.vcf --prefix existing_db

optional arguments:
    -h, --help                      show this help message and exit

    --files [FILES [FILES ...]]      create a db using the specified vcf files (cannot be
                                    used with --folder)

    --folder FOLDER                 create a db using all the vcf files in the folders

    --prefix PREFIX                 the prefix of the output file, default = SVDB

    --upgrade                       create the INS sequence/length table in an existing
                                    database (safe to run on any database; exits with INFO
                                    if already up to date). Optionally combine with --files
                                    or --folder to backfill insertion data from the original
                                    VCFs without rebuilding the full database.

    --pass_only                     only include variants with PASS or . in the FILTER field
                                    (--passonly is a deprecated alias; emits a warning)

    --debug                         enable debug logging to stderr
```

## Export

This module is used to export the variants of the SVDB sqlite database.

Export uses a two-pass clustering approach: DBSCAN-inspired spatial grouping (first pass, controlled by
`--epsilon`/`--min_pts`) followed by overlap/SVLEN/sequence refinement and representative selection
(`--cluster_method`). `--coarse` skips the second pass. See [docs/algorithms.md](docs/algorithms.md)
for a detailed description of all three algorithms.

When the database was built with insertion sequence data (i.e. the INS table is present), insertions are exported with the actual insertion sequence in the ALT column instead of the symbolic `<INS>` allele. For clusters containing multiple samples, the most common sequence across the cluster is used as the representative ALT allele. If the INS table is absent (older database), a warning is emitted and insertions are exported as `<INS>`; run `svdb --build --upgrade --files <original_vcfs> --prefix <existing_db>` to create the INS table and backfill insertion data from the original VCFs.

```text
print a help message
    svdb  --export --help
Export the variants of the database database.db:
    svdb --export --db database.db

optional arguments:
    --no_merge                  skip the merging of variants, print all variants in the db to a vcf file

    --bnd_distance BND_DISTANCE the maximum distance between two similar precise breakpoints (default = 2500)

    --ins_distance INS_DISTANCE the maximum distance to cluster two insertions
                                (default: 25; profile cohort/position_only: 50)

    --ins_svlen_ratio RATIO      minimum SVLEN ratio (min/max) for insertion clustering
                                (default: 0.90; profile cohort: 0.80); requires INS table

    --ins_seq_similarity THRESHOLD
                                minimum Levenshtein sequence similarity (0-1) for insertion clustering;
                                explicit value overrides --data_profile for this parameter
                                (effective default: 0.75); requires INS table

    --data_profile {sample,cohort,position_only}
                                insertion matching preset; requires INS table:
                                  sample:        strict (dist=25, ratio=0.90, sim=0.85) - same individual / technology
                                  cohort:        permissive (dist=50, ratio=0.80, sim=0.75) - cross-individual or cross-caller
                                  position_only: no sequence gate, cluster on position+SVLEN only (dist=50, ratio=0.90)
                                individual --ins_* flags override profile values

    --overlap OVERLAP           the overlap required to merge two events (0 means anything that
                                touches will be merged, 1 means that two events must be identical
                                to be merged), default = 0.8

    --coarse                    skip the second-pass refinement (overlap/SVLEN/sequence gates) and
                                use centroid-based representative selection directly from the DBSCAN
                                spatial clusters. Produces fewer, coarser clusters. Combine with
                                --epsilon/--min_pts to tune the DBSCAN grouping.
                                (--DBSCAN is a deprecated alias; emits a warning)

    --epsilon EPSILON           used together with --coarse; spatial grouping radius in bp
                                (DBSCAN-style epsilon): variants within this distance of each
                                other are candidates for the same cluster (default = 500)

    --min_pts MIN_PTS           used together with --coarse; DBSCAN-style min_pts: minimum number
                                of consecutively-positioned variants that must all fall within
                                --epsilon to seed a cluster. Isolated variants become singletons.
                                (default = 2, meaning any pair within --epsilon forms a cluster)

    --prefix PREFIX             the prefix of the output file, default = same as input

    --memory                    load the database into memory: increases the memory requirements,
                                but lowers the time consumption

    --strip_chr                 strip the 'chr' prefix from chromosome names in the output VCF
                                (e.g. 'chr1' -> '1')

    --samples {on,off}          include per-sample genotype columns (default: on).
                                Use 'off' for sites-only output (analogous to gnomAD --sites-only):
                                the FORMAT and GT columns are omitted; OCC and FRQ are still written

    --max_ins_seq_len N         cap the insertion sequence length used for similarity comparison.
                                Sequences longer than N bp are treated as unknown for clustering
                                purposes (position + SVLEN only); the ALT field falls back to
                                <INS> with SVLEN for those variants. Useful for controlling
                                memory and runtime when the database contains very long sequences.

    --cluster_method {star,union_find}
                                second-pass clustering algorithm applied within each DBSCAN spatial
                                group (default: star).

                                star        greedy star clustering: the highest-degree variant
                                            becomes the representative and claims its direct
                                            neighbours. No transitivity: if A overlaps B and B
                                            overlaps C but A does not overlap C, A and C end up
                                            in separate clusters. Faster; produces more, smaller
                                            clusters.

                                union_find  transitive closure via Union-Find: A-B and B-C are
                                            always merged even if A and C do not overlap directly.
                                            Produces fewer, larger clusters with higher OCC counts.
                                            Prefer this when you want maximal merging across a
                                            cohort; use star when you want conservative merging.

    --workers N                 number of parallel worker processes for the clustering step
                                (default: 0 = use all logical CPUs; 1 = serial).
                                To find your optimal N: time with --workers 1, then increase
                                until wall-clock time stops improving - that plateau is the
                                serial floor (DB fetch, DBSCAN-style grouping, file I/O). On
                                shared systems set N explicitly to avoid monopolising cores.
```

## Query

The query module is used to query one or more structural variant databases. Typically a database is constructed using the build module. However, since this module utilizes the genotype field of the structural variant database vcf to compute the frequency of structural variants, a wide range of files could be used as database. The query module requires a query vcf, as well as a database file (either multisample vcf or SVDB sqlite database):

```text
print a help message
   svdb --query --help
Query a structural variant database, using a vcf file as query:

    svdb --query --query_vcf patient1.vcf --db control_db.vcf

Query multiple databases, using a vcf file as query:

    svdb --query --query_vcf patient1.vcf --db control_db1.vcf,control_db2.vcf --prefix test --in_occ default,Obs --in_frq FRQ,default --out_frq db1_AF,db2_Frq --out_occ db1_AC,db2_Obs

optional arguments:

    -h, --help              show this help message and exit
    --db DB                 path to a db vcf, or a comma separated list of vcfs
    --sqdb SQDB             path to a SVDB sqlite db, or a comma separated list of dbs
    --bedpedb BEDPEDB       path to a SV database of the following format chrA-posA-chrB-posB-type-count-frequency,
                            or a comma separated list of files
    --in_occ IN_OCC         The allele count tag, if used, this tag must be present in the INFO column of the
                            input DB (usually set to AN or OCC). Required if multiple databases are queried.
    --in_frq IN_FRQ         The frequency count tag, if used, this tag must be present in the INFO column of the
                            input DB (usually set to AF or FRQ). Required if multiple databases are queried.
    --out_occ OUT_OCC       the allele count tag, as annotated by SVDB variant (default=OCC).
                            Required if multiple databases are queried.
    --out_frq OUT_FRQ       the tag used to describe the frequency of the variant (default=FRQ).
                            Required if multiple databases are queried.
    --prefix PREFIX         the prefix of the output file, default = print to stdout.
                            Required if multiple databases are queried.
    --bnd_distance BND_DISTANCE
                            the maximum distance between two similar breakpoints (default = 10000)
    --overlap OVERLAP       the overlap required to merge two events (0 means anything that
                            touches will be merged, 1 means that two events must be identical
                            to be merged), default = 0.6
    --ins_distance INS_DISTANCE
                            maximum distance to match two insertions
                            (default: 25; profile cohort/position_only: 50)
                            Applied with --db; also applied with --sqdb when the database
                            contains the INS table; no effect with --bedpedb
    --ins_svlen_ratio INS_SVLEN_RATIO
                            minimum SVLEN ratio (min/max) required to match two insertions
                            with known length (default: 0.90; profile cohort: 0.80)
                            Applied with --db; also applied with --sqdb when the database
                            contains the INS table; no effect with --bedpedb
    --ins_seq_similarity THRESHOLD
                            minimum Levenshtein sequence similarity (0-1) required to match
                            two insertions with known sequence; explicit value overrides
                            --data_profile for this parameter (effective default: 0.75)
                            Applied with --db; also applied with --sqdb when the database
                            contains the INS table; no effect with --bedpedb
    --data_profile {sample,cohort,position_only}
                            insertion matching preset:
                              sample:        strict (dist=25, ratio=0.90, sim=0.85) - same individual / technology
                              cohort:        permissive (dist=50, ratio=0.80, sim=0.75) - cross-individual or cross-caller
                              position_only: no sequence gate, match on position+SVLEN only (dist=50, ratio=0.90)
                            individual --ins_* flags override profile values
                            Applied with --db; also applied with --sqdb when the database
                            contains the INS table; no effect with --bedpedb

    --max_ins_seq_len N     sequences longer than N bp are excluded from sequence similarity
                            matching and fall back to position+SVLEN; off by default.
                            500 or 1000 is recommended for large databases to reduce runtime.

    --memory                load the database into memory: increases the memory requirements,
                            but lowers the time consumption (may only be used with sqdb)
    --no_var                count overlapping variants of different type as hits in the db
```

## Merge

The merge module merges variants within one or more vcf files. This could be used to either merge the output of multiple callers, or to merge variants that are called multiple times due to noise or some other error:

```text
print a help message:
   python SVDB.py --merge --help
merge vcf files:
    svdb --merge --vcf patient1_lumpy.vcf patient1_cnvnator.vcf patient1_TIDDIT.vcf > patient1_merged_callers.vcf

Similar variants will be merged, and presented according to the order of the input vcf files. I.e If lumpy and cnvnator calls the same variant in the top example,
the variant will be printed as the lumpy call. In most cases, the order should be set according to the accuracy or detail of the info field of the different callers.
The order could also be set using the --priority flag:
    svdb --merge --vcf patient1_lumpy.vcf:one patient1_cnvnator.vcf:2 patient1_TIDDIT.vcf:tiddit --priority tiddit,2,one > patient1_merged_callers.vcf

In this example, tiddit will have the highest order, cnvnator second etc.

optional arguments:
    -h, --help                      show this help message and exit

    --bnd_distance BND_DISTANCE     the maximum distance between two similar precise breakpoints
                                    (default = 2000)

    --overlap OVERLAP               the overlap required to merge two events (0 means
                                    anything that touches will be merged, 1 means that two
                                    events must be identical to be merged), default = 0.95

    --ins_distance INS_DISTANCE     maximum distance to merge two insertions
                                    (default: 25; profile cohort/position_only: 50)

    --ins_svlen_ratio INS_SVLEN_RATIO
                                    minimum SVLEN ratio (min/max) required to merge two
                                    insertions with known length (default: 0.90; profile cohort: 0.80)

    --ins_seq_similarity THRESHOLD  minimum Levenshtein sequence similarity (0-1) required to
                                    merge two insertions with known sequence; explicit value
                                    overrides --data_profile for this parameter
                                    (effective default: 0.75)

    --data_profile {sample,cohort,position_only}
                                    insertion matching preset:
                                      sample:        strict (dist=25, ratio=0.90, sim=0.85) - same individual / technology
                                      cohort:        permissive (dist=50, ratio=0.80, sim=0.75) - cross-individual or cross-caller
                                      position_only: no sequence gate, merge on position+SVLEN only (dist=50, ratio=0.90)
                                    individual --ins_* flags override profile values

    --priority                      prioritise the input vcf files

    --no_intra                      no merging of variants within the same vcf

    --no_var                        variants of different type will be merged

    --pass_only                     merge only variants labeled PASS

    --same_order                    assume that the samples are ordered the same way (skip
                                    reordering and merging of the sample columns)
```

## For developers

Runtime dependencies are pinned via [pip-tools](https://pip-tools.readthedocs.io/). Edit `requirements.in`, then regenerate:

```bash
pip-compile requirements.in --output-file requirements.txt --strip-extras
```

Dev tools (`requirements-dev.txt`) are intentionally unpinned — they span Python 3.9–3.14 where transitive pins would differ per version.

Install development dependencies:

```bash
pip install -r requirements-dev.txt
```

Run tests (includes ruff linting and mypy type checking):

```bash
pytest
```

Run ruff or mypy standalone:

```bash
ruff check svdb/
mypy svdb/ --ignore-missing-imports
```

Configuration lives in `pyproject.toml` (build system, ruff, pytest settings). The legacy `setup.py` is retained only for optional Cython compilation of `merge_vcf_module_cython`.

See [docs/architecture.md](docs/architecture.md) for a module overview and data flow diagrams.

## Profiling

A cProfile-based profiling harness lives in `scripts/profile_svdb.py`. It runs a standard battery of commands (merge, build, export, query) on real VCF data and prints per-function timing.

Set up a local config file (gitignored):

```bash
cp scripts/profile_config.toml.example scripts/profile_config.toml
# fill in your VCF paths and caller names
```

Then run:

```bash
python scripts/profile_svdb.py               # default: top 15 functions, sorted by cumulative time
python scripts/profile_svdb.py --top 20 --sort tottime
python scripts/profile_svdb.py --config /path/to/my_config.toml
```

The script always profiles the local checkout (not any installed package), so it is safe to use during optimisation work.
