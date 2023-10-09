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

external databases are run like this:

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
	
	git clone https://github.com/J35P312/SVDB.git
	cd SVDB
	pip install -e .

SVDB is available on singularity:

	singularity pull shub://J35P312/SVDB

# Modules:
SVDB consists of modules that are used to build, query, annotate, export, and analyse structural variant databases. These are the modules:

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
This module is used to export the variants of the SVDB sqlite database. The variants of the sqlite svdb database is clustered using one out of three algorithms, overlap or DBSCAN.
 
    print a help message
        svdb  --export --help  
    Export the variants of the database database.db:
        svdb --export --db database.db

    optional arguments:
        --no_merge            skip the merging of variants, print all variants in the db to a vcf file

         --bnd_distance BND_DISTANCE  the maximum distance between two similar precise breakpoints(default = 2500)
 
         --overlap OVERLAP     the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.8

         --DBSCAN              use dbscan to cluster the variants, overides the overlap based clustering algorithm

         --epsilon EPSILON     used together with --DBSCAN; sets the epsilon parameter(default = 500bp)

         --min_pts MIN_PTS     the min_pts parameter(default = 2

         --prefix PREFIX       the prefix of the output file, default = same as input

          --memory              load the database into memory: increases the memory requirements, but lowers the time consumption

## Query
The query module is used to query one or more structural variant databases. Typically a database is constructed using the build module. However, since this module utilize the genotype field of the structural variant database vcf to compute the frequency of structural variants, a wide range of files could be used as database. The query module requires a query vcf, as well as a database file(either multisample vcf or SVDB sqlite database):

    print a help message
       svdb --query --help
    Query a structural variant database, using a vcf file as query:

        svdb --query --query_vcf patient1.vcf --db control_db.vcf

    Query multiple databases, using a vcf file as query:
    
        svdb --query --query_vcf patient1.vcf --db control_db1.vcf,control_db2.vcf --prefix test --in_occ default,Obs --in_frq FRQ,default --out_frq db1_AF,db2_Frq --out_occ db1_AC,db2_Obs

    optional arguments:

	-h, --help            show this help message and exit

	--db DB				path to a db vcf, or a comma separated list of vcfs
	--sqdb SQDB			path to a SVDB sqlite db, or a comma separated list of dbs
	--bedpedb BEDPEDB		path to a SV database of the following format chrA-posA-chrB-posB-type-count-frequency, or a comma separated list of files
	--in_occ IN_OCC			The allele count tag, if used, this tag must be present in the INFO column of the input DB(usually set to AN or OCC). This parameter is required if multiple databases are queried. 
	--in_frq IN_FRQ			The frequency count tag, if used, this tag must be present in the INFO column of the input DB(usually set to AF or FRQ). This parameter is required if multiple databases are queried.
	--out_occ OUT_OCC		the allele count tag, as annotated by SVDB variant(default=OCC). This parameter is required if multiple databases are queried.
	--out_frq OUT_FRQ		the tag used to describe the frequency of the variant(default=FRQ). This parameter is required if multiple databases are queried.
	--prefix PREFIX			the prefix of the output file, default = print to stdout. Required if multiple databases are queried. 
	--bnd_distance BND_DISTANCE	the maximum distance between two similar breakpoints(default = 10000)
	--overlap OVERLAP		the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6
	--memory 			load the database into memory: increases the memory requirements, but lowers the time consumption(may only be used with sqdb)
	--no_var			count overlapping variants of different type as hits in the db

    
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

## Annotate
The annotate module works similarily to the query module. However, the annotate module is used to annotate a query vcf with any selected INFO tags from a vcf database. The annotate module does not support db or bedbpe files.

    print a help message:
       python SVDB.py --annotate --help
    merge vcf files:
        svdb --annotate --query_vcf input.vcf --db frequency_db.vcf --in_tags AF SVTYPE --out_tag LOCALAF TYPE_of_SV > annotated.vcf

In this example, AF and SVTYPE are extracted from frequency_db.vcf, the tags will be stored as LOCALAF and TYPE_of_SV in the info column of matching variants in the annotated.vcf file.

optional arguments:
  -h, --help            show this help message and exit
  --annotate            query a db
  --query_vcf QUERY_VCF
                        a vcf used to query the db
  --db DB               path to a database vcf
  --in_tags [IN_TAGS [IN_TAGS ...]]
                        list of tags to extract from the INFO field of the vcf
                        database
  --out_tags [OUT_TAGS [OUT_TAGS ...]]
                        the name of the tags as written to the query vcf
                        Default=same as in_tags
  --prefix PREFIX       the prefix of the output file, default = print to
                        stdout. Required, if multiple databases are queried
  --bnd_distance BND_DISTANCE
                        the maximum distance between two similar
                        breakpoints(default = 10000)
  --ins_distance INS_DISTANCE
                        the maximum distance to merge two insertions(default =
                        50)
  --overlap OVERLAP     the overlap required to merge two events(0 means
                        anything that touches will be merged, 1 means that two
                        events must be identical to be merged), default = 0.6
  --no_var              count overlaping variants of different type as hits in
 
