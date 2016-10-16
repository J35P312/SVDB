from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(["SVDB_build_module.py","SVDB_overlap_module.py","readVCF.py","SVDB_merge_vcf_module_cython.py","SVDB_hist_module.py","SVDB_query_module_SQLITE.py"])
)
