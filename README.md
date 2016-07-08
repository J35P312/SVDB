# SVDB
SVDB is a toolkit for constructing and querying structural variant databases. The databases are constructed using the output vcf files from structural variant callers such as TIDDIT, Manta, Fermikit or Delly.
The thousand genomes structural variant calls may also be used as a database:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/

#Install:
Dependencies: SVDB has been tested on python 2.7.11, only standard modules are required.
Optional: the code may be compiled using cython, this will speed up the database construction(this requires cython):

python setup.py build_ext --inplace

#modules:
SVDB consists of five separate modules that are used to manage, query and create structural variant databases. These are the modules:

    Build:
      This module is used to construct structural variant databases from vcf files. It is activated using the folowing command:
        python svdb.py --build
      for more information type:
        python SVDB.py --build --help
        
    Hist:
      This module is used to compare structural variant vcf files, either by generating a similarity matrix, or by creating histograms of the efficency of databases of different sizes(based on input vcf files). The module is activated through this command:
        python svdb.py --hist
      for more information type:
        python SVDB.py --hist --help
    Query:
    The query module is used to query a structural variant database. Typically a database is constructed using the build module. However, since this module utilize the genotype field of the sructural variant database vcf to compute the frequency of structural variants, a wide range of files could be used as database. The query module requires a query vcf, as well as a database vcf.
      python SVDB.py --query
    for more info, type
       python SVDB.py --query --help
    Purge:
      The purge module is used to remove entries from a database. These entries could be sensitive information such as disease causing variants. THe purge module is run using this command:
        python SVDB.py --purge
      for more info, type
       python SVDB.py --purge --help
       
    Merge:
      The merge module merges variants within one or more vcf files. This could be used to either merge the output of multiple callers, or to merge variants that are called multiple times due to noise or some other error.
        python SVDB.py --merge
      for more info, type
       python SVDB.py --merge --help
