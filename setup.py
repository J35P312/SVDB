from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(["SVDB_overlap_module.py","SVDB_build_module_cython.py","readVCF.py","SVDB_merge_vcf_module_cython.py"])
)
