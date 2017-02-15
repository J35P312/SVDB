# SVDB
SVDB is a toolkit for constructing and querying structural variant databases. The databases are constructed using the output vcf files from structural variant callers such as TIDDIT, Manta, Fermikit or Delly.
The thousand genomes structural variant calls may also be used as a database:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/

#Install:
Dependencies: SVDB has been tested on python 2.7.11, and requires sciKit-learn as well as numpy.
These two packages may be installed using pip:

pip install numpy

pip install -U scikit-learn

Optional: the code may be compiled using cython, this will speed up the database construction(this requires cython):

python setup.py install

#modules:
SVDB consists of five separate modules that are used to manage, query and create structural variant databases. These are the modules:

Build: This module is used to construct structural variant databases from vcf files. The database may then be queried to compute the frequency of structural variants. These are the commands used to construct a structural variation database:
    
    print a help message
        svdb  --build --help  
    Construct a database, from a set of vcf files:
        svdb --build --vcf sample1.vcf sample2.vcf sample3.vcf
    construct a database from vcf files stored in a folder
        svdb --build --folder SV_analysis_folder/
        
    optional arguments:
        -h, --help                      show this help message and exit

        --no_merge                      skip the clustering of variants
  
        --ci                            overides overlap and bnd_distance,determine hits based
                                        on the confidence interval of the position fo the
                                        variants(0 if no CIPOS or CIEND is vailable)
                                        
        --bnd_distance BND_DISTANCE     the maximum distance between two similar precise
                                        breakpoints(default = 2500) 
                        
        --overlap OVERLAP               the overlap required to merge two events(0 means
                                        anything that touches will be merged, 1 means that two
                                        events must be identical to be merged), default = 0.8
                                        
        --files [FILES [FILES ...]]      create a db using the specified vcf files(cannot be
                                         used with --folder)
                        
        --folder FOLDER                  create a db using all the vcf files in the folders
        
        --prefix PREFIX                  the prefix of the output file, default = SVDB

        
Hist: This module is used to compare structural variant vcf files, either by generating a similarity matrix, or by creating histograms of the efficency of databases of different sizes(based on input vcf files):

    print a help message
        python SVDB.py --hist --help
    Create histograms of different sizes, and compute their efficiency:
        python --hist --sample_hist -folder input
    Create a similarity matrix of the selected sampes:
        python --hist --similarity_matrix -folder input
    
    optional arguments:
    
        -h, --help                      show this help message and exit
        
        --files [FILES [FILES ...]]     input vcf files(cannot be used with folder)
         
        --k [K [K ...]]                 the sizes of the sampled databases
                                        default = n=10*i < samples(used with sample_hist)
        
        --n N                           the number of iterations,default=100(used with sample_hist)
  
        --bnd_distance BND_DISTANCE     the maximum distance between two similar precise
                                        breakpoints(default = 10000)
                                        
        --overlap OVERLAP               the overlap required to merge two events(0 means
                                        anything that touches will be merged, 1 means that two
                                        events must be identical to be merged), default = 0.6
        
        --ci                            overides overlap and bnd_distance,determine hits based
                                        on the confidence interval of the position of the
                                        variants(0 if no CIPOS or CIEND is vailable)

Query: The query module is used to query a structural variant database. Typically a database is constructed using the build module. However, since this module utilize the genotype field of the sructural variant database vcf to compute the frequency of structural variants, a wide range of files could be used as database. The query module requires a query vcf, as well as a database vcf:

    print a help message
       python SVDB.py --query --help
    Query a structural variant database, using a vcf file as query:  
        svdb --query --query_vcf patient1.vcf --db control_db.vcf
    
    optional arguments:
        -h, --help                      show this help message and exit

        --hit_tag HIT_TAG               the tag used to describe the number of hits within the
                                        info field of the output vcf(default=OCC)
                                        
        --frequency_tag FREQUENCY_TAG   the tag used to describe the frequency of the
                                        variant(defualt=FRQ)
                        
        --prefix PREFIX                 the prefix of the output file, default = print to stdout

        --bnd_distance BND_DISTANCE     the maximum distance between two similar precise breakpoints
                                        (default = 10000)
                        
                        
        --overlap OVERLAP               the overlap required to merge two events(0 means
                                        anything that touches will be merged, 1 means that two
                                        events must be identical to be merged), default = 0.6
                                        
        --no_var                        count overlaping variants of different type as hits
        
        --invert                        invert the sorting order so that high frequency
                                        samples are present on top of the output vcf
                              
        --ci                            overides overlap and bnd_distance,determine hits based
                                        on the confidence interval of the position fo the
                                        variants(0 if no CIPOS or CIEND is vailable)

Purge: The purge module is used to remove entries from a database:

    print a help message:
       python SVDB.py --purge --help
    Delete a sample from a DB, the sample id should be the same as the id written in the format columns of the db:
        svdb --purge --sample patient2 --db my_svdb.vcf > cleaned_db.vcf
    Delete variants from a DB, the variants should be stored in a standard structural variant format:
        svdb --purge --vcf delete_these_variants.vcf --db my_svdb.vcf > cleaned_db.vcf
    
    optional arguments:
        -h, --help                      show this help message and exit
        
        --bnd_distance BND_DISTANCE     the maximum distance between two similar precise breakpoints
                                        (default = 10000)
                        
        --overlap OVERLAP               the overlap required to merge two events(0 means
                                        anything that touches will be merged, 1 means that two
                                        events must be identical to be merged), default = 0.6
                              
        --ci                            overides overlap and bnd_distance,determine hits based
                                        on the confidence interval of the position fo the
                                        variants(0 if no CIPOS or CIEND is vailable)

    
    
Merge: The merge module merges variants within one or more vcf files. This could be used to either merge the output of multiple callers, or to merge variants that are called multiple times due to noise or some other error:

    print a help message:
       python SVDB.py --merge --help
    merge vcf files:
        svdb --merge --vcf patient1_lumpy.vcf patient1_cnvnator.vcf patient1_TIDDIT.vcf > patient1_merged_callers.vcf 

    optional arguments:
        -h, --help                      show this help message and exit
        
        --bnd_distance BND_DISTANCE     the maximum distance between two similar precise breakpoints
                                        (default = 10000)
                        
        --overlap OVERLAP               the overlap required to merge two events(0 means
                                        anything that touches will be merged, 1 means that two
                                        events must be identical to be merged), default = 0.6
                              
        --ci                            overides overlap and bnd_distance,determine hits based
                                        on the confidence interval of the position fo the
                                        variants(0 if no CIPOS or CIEND is vailable)
                                        
        --no_intra                       no merging of variants within the same vcf
        
        --no_var                        variants of different type will be merged
        
        --pass_only                     merge only variants labeled PASS
