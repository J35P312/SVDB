# SVDB
SVDB is a toolkit for constructing and querying structural variant databases. The databases are constructed using the output vcf files from structural variant callers such as TIDDIT, Manta, Fermikit or Delly.
The thousand genomes structural variant calls may also be used as a database:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/

#Install:
the code may be compiled using cython, this will speed up the database construction(this requires cython):

python setup.py build_ext --inplace

#Commands
for more info, run the SVDB script:

python SVDB.py --help

for info on a specific module type the module name, example:

python SVDB.py --build --help
