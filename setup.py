from setuptools import setup

try:
    from Cython.Build import cythonize
    has_cython = True
except ImportError:
    has_cython = False

if has_cython:
    ext_modules = cythonize([
        "SVDB_build_module.py",
        "SVDB_overlap_module.py",
        "SVDB_readVCF.py",
        "SVDB_merge_vcf_module_cython.py",
        "SVDB_hist_module.py",
        "SVDB_query_module.py"])
else:
    ext_modules = []

setup(
    name = 'svdb',
    version = '0.1.0',
    ext_modules = ext_modules,
    packages = ['svdb'],
    install_requires = ['numpy', 'scikit-learn', 'scipy'],
    entry_points = {'console_scripts': ['svdb = svdb.__main__:main']},
)
