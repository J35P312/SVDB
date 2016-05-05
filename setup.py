from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(["SVDB_overlap_module.pyx","SVDB_build_module_cython.pyx","readVCF.pyx","SVDB_merge_vcf_module_cython.pyx"])
)
