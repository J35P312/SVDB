# SVDB
Dependencies:cython

Install:

the code may be compiled using cython, this will speed up the database construction:

python setup.py build_ext --inplace

for more info, run the SVDB script:

python SVDB.py --help

for info on a specific module type the module name, example:

python SVDB.py --build --help
