# SVDB
SVDB is a toolkit for constructing and querying structural variant databases. The databases are constructed using the output vcf files from structural variant callers such as TIDDIT, Manta, Fermikit or Delly.
SVDB may also be used to merge SV  vcf files from multiple callers or individuals.


# Supported public databases
SVDB query supports public databases such as thousand genomes SV map and Gnomad SV, as well as most multisample SV vcf files

The thousand genomes SV database:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/

The swegen SVDB:

https://swefreq.nbis.se/

The GNOMAD SV database:

https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2_sv.sites.vcf.gz

external databses are run like this:

```bash
svdb --query \
     --query_vcf /home/jesper/vcf/6_pair_limit/P2109_120.clean.dedup.recal_FindSV.vcf \
     --out_occ GNOMAD_AC \
     --out_frq GNOMAD_AF \
     --in_occ AN \
     --out_frq AF \
     --db /home/jesper/Downloads/gnomad_sv/gnomad_v2_sv.sites.vcf
```

here the `AF` and `AN` are the allele frequency tags of the database, the `AF` is a float, and `AN` is an integer. These tags will be added to the annotated output vcf, and named `GNOMAD_AC`, `GNOMAD_AF`.

# Install:
Dependencies: SVDB has been tested on python 2.7.11 and python 3.6, and requires numpy.
SVDB is installed using the following command

	pip install -e .

SVDB is available on singularity:

	singularity pull shub://J35P312/SVDB

# Modules:
SVDB consists of modules that are used to build, query, export, and analyse structural variant databases. These are the modules:

## Build
This module is used to construct structural variant databases from vcf files. The database may then be queried to compute the frequency of structural variants, or exported into a vcf file. These are the commands used to construct a structural variation database:
    
    print a help message
        svdb  --build --help  
    Construct a database, from a set of vcf files:
        svdb --build --vcf sample1.vcf sample2.vcf sample3.vcf
    construct a database from vcf files stored in a folder
        svdb --build --folder SV_analysis_folder/
        
    optional arguments:
        -h, --help                      show this help message and exit

        --files [FILES [FILES ...]]      create a db using the specified vcf files(cannot be
                                         used with --folder)
                        
        --folder FOLDER                  create a db using all the vcf files in the folders
        
        --prefix PREFIX                  the prefix of the output file, default = SVDB


## Export
This module is used to export the variants of the SVDB sqlite database. The variants of the sqlite svdb database is clustered using one out of three algorihms, overlap or DBSCAN.
 
    print a help message
        svdb  --export --help  
    Export the variants of the database database.db:
        svdb --export --db database.db

    optional arguments:
        --no_merge            skip the merging of variants, print all variants in the db to a vcf file

         --bnd_distance BND_DISTANCE  the maximum distance between two similar precise breakpoints(default = 2500)
 
         --overlap OVERLAP     the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.8

         --DBSCAN              use dbscan to cluster the variants, overides the overlap based clustering algoritm

         --epsilon EPSILON     used together with --DBSCAN; sets the epsilon paramter(default = 500bp)

         --min_pts MIN_PTS     the min_pts parameter(default = 2

         --prefix PREFIX       the prefix of the output file, default = same as input

          --memory              load the database into memory: increases the memory requirements, but lowers the time consumption

## Query
The query module is used to query a structural variant database. Typically a database is constructed using the build module. However, since this module utilize the genotype field of the sructural variant database vcf to compute the frequency of structural variants, a wide range of files could be used as database. The query module requires a query vcf, as well as a database file(either multisample vcf or SVDB sqlite database):

    print a help message
       svdb --query --help
    Query a structural variant database, using a vcf file as query:

        svdb --query --query_vcf patient1.vcf --db control_db.vcf

    optional arguments:

	-h, --help            show this help message and exit

	--db DB				path to a db vcf
	--sqdb SQDB			path to a SVDB sqlite db
	--bedpedb BEDPEDB		path to a SV database of the following format chrA-posA-chrB-posB-type-count-frequency
	--in_occ IN_OCC			The allele count tag, if used, this tag must be present in the INFO column of the input DB(usually set to AN or OCC)
	--in_frq IN_FRQ			The frequency count tag, if used, this tag must be present in the INFO column of the input DB(usually set to AF or FRQ)
	--out_occ OUT_OCC		the allle count tag, as annotated by SVDB variant(defualt=OCC)
	--out_frq OUT_FRQ		the tag used to describe the frequency of the variant(defualt=FRQ)
	--prefix PREFIX			the prefix of the output file, default = print to stdout
	--bnd_distance BND_DISTANCE	the maximum distance between two similar breakpoints(default = 10000)
	--overlap OVERLAP		the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6
	--memory 			load the database into memory: increases the memory requirements, but lowers the time consumption(may only be used with sqdb)
	--no_var			count overlaping variants of different type as hits in the db

    
## Merge
The merge module merges variants within one or more vcf files. This could be used to either merge the output of multiple callers, or to merge variants that are called multiple times due to noise or some other error:

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
                                        (default = 10000)
                        
        --overlap OVERLAP               the overlap required to merge two events(0 means
                                        anything that touches will be merged, 1 means that two
                                        events must be identical to be merged), default = 0.6

        --priority                      prioritise the input vcf files
                                                                      
        --no_intra                      no merging of variants within the same vcf
        
        --no_var                        variants of different type will be merged
        
        --pass_only                     merge only variants labeled PASS

	--same_order			assume that the samples are ordered the same way (skip reordering and merging of the sample columns). 
