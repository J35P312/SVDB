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

SVDB is available as a container on [BioContainers](https://quay.io/repository/biocontainers/svdb?tab=tags&tag=).

## Global options

These flags apply to all subcommands and must be placed before the subcommand name:

```text
    --debug     enable debug logging to stderr
                e.g. svdb --debug --build --files sample.vcf
```

## Modules

SVDB consists of modules that are used to build, query, export, and analyse structural variant databases. These are the modules:

## Build

This module is used to construct structural variant databases from vcf files. The database may then be queried to compute the frequency of structural variants, or exported into a vcf file. These are the commands used to construct a structural variation database:

Sample names are taken from the VCF header's sample column (the named `FORMAT` columns).
For VCFs with no sample columns (INFO-only format), the filename stem is used instead.

```text
    svdb --build --help
    svdb --build --files sample1.vcf sample2.vcf sample3.vcf
    svdb --build --folder SV_analysis_folder/
    svdb --build --upgrade --files sample1.vcf sample2.vcf --prefix existing_db   # add INS table + backfill insertion sequences from original VCFs

  input (one of required):
    --files [FILES ...]             vcf files to build the db from (cannot be used with --folder)
    --folder FOLDER                 use all vcf files in the given folder

  output:
    --prefix PREFIX                 prefix for the output file (default: SVDB)
    --pass_only                     only include variants with PASS or . in the FILTER field
                                    (--passonly is a deprecated alias; emits a warning)

  upgrade existing db:
    --upgrade                       create the INS table and backfill insertion sequences from the
                                    provided VCFs; schema-only, no existing SVDB rows are changed.
                                    Exits with INFO if the INS table already exists.
                                    Warns (WARNING) for DB samples with no matching VCF provided;
                                    logs (INFO) VCF samples that have no entries in the database.
                                    Requires --files or --folder.

  storage:
    --max_ins_seq_len N             cap on stored insertion sequence length (bp); sequences
                                    longer than N are stored with NULL sequence but retain
                                    SVLEN for length-ratio matching (default: no limit)
```

## Export

This module is used to export the variants of the SVDB sqlite database.

Export uses a two-pass clustering approach: DBSCAN-inspired spatial grouping (first pass, controlled by
`--epsilon`/`--min_pts`) followed by overlap/SVLEN/sequence refinement and representative selection
(`--cluster_method`). `--coarse` skips the second pass. See [docs/algorithms.md](docs/algorithms.md)
for a detailed description of all three algorithms.

When the database was built with insertion sequence data (i.e. the INS table is present), insertions are exported with the actual insertion sequence in the ALT column instead of the symbolic `<INS>` allele. For clusters containing multiple samples, the most common sequence across the cluster is used as the representative ALT allele. If any cluster member lacks a sequence — because its insertion was longer than `--max_ins_seq_len` at build time, or because it was called as a symbolic `<INS>` allele — the entire cluster is exported as `<INS>` with SVLEN taken from the most common stored length across members; a mixed cluster cannot be faithfully represented by a single sequence. If the INS table is absent (older database), a warning is emitted and insertions are exported as `<INS>`; run `svdb --build --upgrade --files <original_vcfs> --prefix <existing_db>` to create the INS table and backfill insertion data from the original VCFs.

```text
    svdb --export --help
    svdb --export --db database.db

  input (required):
    --db DB                     the SQLite database to export

  output:
    --prefix PREFIX             prefix for the output file (default: same as input)
    --no_merge                  skip merging; print all variants as-is
    --strip_chr                 strip the 'chr' prefix from chromosome names in the output
                                (e.g. 'chr1' -> '1'); names are stored as-is in the db
    --samples {on,off}          include per-sample genotype columns (default: on);
                                use 'off' for sites-only output (FORMAT/GT omitted, OCC/FRQ kept)

  algorithm — SV matching:
    --bnd_distance BND_DISTANCE maximum distance between two similar precise breakpoints (default: 2500)
    --overlap OVERLAP           minimum reciprocal overlap to merge two events;
                                must be in [0.0, 1.0] (0 = anything touching; 1 = identical) (default: 0.8)

  algorithm — insertion matching (requires INS table):
    --data_profile {sample,cohort,position_only}
                                preset for all insertion parameters; individual --ins_* flags override:
                                  sample:        strict  (dist=25, ratio=0.90, sim=0.85) - same individual / technology
                                  cohort:        permissive (dist=50, ratio=0.80, sim=0.75) - cross-individual or cross-caller
                                  position_only: no sequence gate (dist=50, ratio=0.90)
    --ins_distance INS_DISTANCE maximum distance to cluster two insertions
                                (default: 25; profile cohort/position_only: 50)
    --ins_svlen_ratio RATIO     minimum SVLEN ratio (min/max) for insertion clustering
                                must be in [0.0, 1.0] (default: 0.90; profile cohort: 0.80)
    --ins_seq_similarity THRESHOLD
                                minimum Levenshtein sequence similarity; must be in [0.0, 1.0];
                                explicit value overrides --data_profile (effective default: 0.75)

  algorithm — clustering:
    --coarse                    skip second-pass refinement; use centroid from DBSCAN spatial
                                clusters directly. Produces fewer, coarser clusters.
                                (--DBSCAN is a deprecated alias; emits a warning)
    --epsilon EPSILON           DBSCAN-style spatial grouping radius in bp (default: 500)
    --min_pts MIN_PTS           DBSCAN-style min_pts: minimum variants within --epsilon to seed
                                a cluster; isolated variants become singletons; must be a whole
                                number ≥ 1 (default: 2, meaning any pair within --epsilon forms a cluster)
    --cluster_method {star,union_find}
                                second-pass clustering algorithm (default: star):
                                  star        greedy; highest-degree variant claims neighbours.
                                              No transitivity. A variant overlapping multiple
                                              representatives appears in each cluster — OCC
                                              reflects all groups it genuinely belongs to.
                                  union_find  transitive closure; A-B + B-C → {A,B,C}.
                                              Exclusive membership; fewer, larger clusters.

  performance:
    --max_ins_seq_len N         sequences longer than N bp fall back to position+SVLEN
                                (default: 1000); use 0 for no cap (all sequences compared)
    --memory                    load the db into memory: higher memory use, lower export time
    --workers N                 parallel worker processes (default: 0 = all logical CPUs; 1 = serial)
```

## Query

The query module is used to query one or more structural variant databases. Typically a database is constructed using the build module. However, since this module utilizes the genotype field of the structural variant database vcf to compute the frequency of structural variants, a wide range of files could be used as database. The query module requires a query vcf, as well as a database file (either multisample vcf or SVDB sqlite database):

```text
    svdb --query --help
    svdb --query --query_vcf patient1.vcf --db control_db.vcf
    svdb --query --query_vcf patient1.vcf --db control_db1.vcf,control_db2.vcf \
         --prefix test --in_occ default,Obs --in_frq FRQ,default \
         --out_frq db1_AF,db2_Frq --out_occ db1_AC,db2_Obs

  input:
    --query_vcf VCF   (required) query vcf file
    --db DB                     db vcf, or a comma-separated list (no effect on --bedpedb)
    --sqdb SQDB                 SVDB sqlite db, or a comma-separated list
    --bedpedb BEDPEDB           SV db in chrA-posA-chrB-posB-type-count-frequency format,
                                or a comma-separated list
                                (at least one of --db / --sqdb / --bedpedb is required)

  output:
    --prefix PREFIX             prefix for the output file (default: print to stdout);
                                required when querying multiple databases
    --out_occ OUT_OCC           output tag for allele count (default: OCC)
    --out_frq OUT_FRQ           output tag for allele frequency (default: FRQ)
    --in_occ IN_OCC             allele count tag in the db INFO column (usually OCC or AN);
                                required when querying multiple databases
    --in_frq IN_FRQ             allele frequency tag in the db INFO column (usually FRQ or AF);
                                required when querying multiple databases

  algorithm — SV matching:
    --bnd_distance BND_DISTANCE maximum distance between two similar breakpoints (default: 10000)
    --overlap OVERLAP           minimum reciprocal overlap to match two events;
                                must be in [0.0, 1.0] (0 = anything touching; 1 = identical) (default: 0.6)
    --max_frq MAX_FRQ           only report variants with frequency at or below this value
                                (default: 1, i.e. all variants)
    --no_var                    count overlapping variants of different type as hits in the db

  algorithm — insertion matching (--db and --sqdb with INS table; no effect with --bedpedb):
    --data_profile {sample,cohort,position_only}
                                preset for all insertion parameters; individual --ins_* flags override:
                                  sample:        strict  (dist=25, ratio=0.90, sim=0.85) - same individual / technology
                                  cohort:        permissive (dist=50, ratio=0.80, sim=0.75) - cross-individual or cross-caller
                                  position_only: no sequence gate (dist=50, ratio=0.90)
    --ins_distance INS_DISTANCE maximum distance to match two insertions
                                (default: 25; profile cohort/position_only: 50)
    --ins_svlen_ratio RATIO     minimum SVLEN ratio (min/max) for insertions with known length
                                must be in [0.0, 1.0] (default: 0.90; profile cohort: 0.80)
    --ins_seq_similarity THRESHOLD
                                minimum Levenshtein sequence similarity; must be in [0.0, 1.0];
                                explicit value overrides --data_profile (effective default: 0.75)

  performance:
    --max_ins_seq_len N         sequences longer than N bp fall back to position+SVLEN
                                (default: 1000); use 0 for no cap (all sequences compared)
    --memory                    load the db into memory: higher memory use, lower query time
                                (sqdb only)
```

## Merge

The merge module merges variants within one or more vcf files. This could be used to either merge the output of multiple callers, or to merge variants that are called multiple times due to noise or some other error:

```text
    svdb --merge --help
    svdb --merge --vcf patient1_lumpy.vcf patient1_cnvnator.vcf patient1_TIDDIT.vcf
    svdb --merge --vcf patient1_lumpy.vcf:one patient1_cnvnator.vcf:2 patient1_TIDDIT.vcf:tiddit \
         --priority tiddit,2,one

Variants are merged and output in the order of the input files (first file takes precedence).
Use --priority to override the order explicitly.

  input (required):
    --vcf VCF [VCF ...]         input vcf files to merge

  input control:
    --priority ORDER            prioritise input files; format: --vcf a.vcf:1 b.vcf:2 --priority 2,1
    --pass_only                 merge only variants labeled PASS
    --no_intra                  skip merging of variants within the same vcf
    --same_order                assume sample columns are in the same order across all input files
    --no_tag                    do not add VARID and set entries to the INFO field
                                (--notag is a deprecated alias; emits a warning)

  algorithm — SV matching:
    --bnd_distance BND_DISTANCE maximum distance between two similar precise breakpoints (default: 2000)
    --overlap OVERLAP           minimum reciprocal overlap to merge two events;
                                must be in [0.0, 1.0] (0 = anything touching; 1 = identical) (default: 0.95)
    --no_var                    variants of different type will be merged

  algorithm — insertion matching:
    --data_profile {sample,cohort,position_only}
                                preset for all insertion parameters; individual --ins_* flags override:
                                  sample:        strict  (dist=25, ratio=0.90, sim=0.85) - same individual / technology
                                  cohort:        permissive (dist=50, ratio=0.80, sim=0.75) - cross-individual or cross-caller
                                  position_only: no sequence gate (dist=50, ratio=0.90)
    --ins_distance INS_DISTANCE maximum distance to merge two insertions
                                (default: 25; profile cohort/position_only: 50)
    --ins_svlen_ratio RATIO     minimum SVLEN ratio (min/max) for insertions with known length
                                must be in [0.0, 1.0] (default: 0.90; profile cohort: 0.80)
    --ins_seq_similarity THRESHOLD
                                minimum Levenshtein sequence similarity; must be in [0.0, 1.0];
                                explicit value overrides --data_profile (effective default: 0.75)

  performance:
    --max_ins_seq_len N         sequences longer than N bp fall back to position+SVLEN
                                (default: 1000); use 0 for no cap (all sequences compared)
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
