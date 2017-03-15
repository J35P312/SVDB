from setuptools import setup

try:
    from Cython.Build import cythonize
    has_cython = True
except ImportError:
    has_cython = False

if has_cython:
    ext_modules = cythonize([
        "svdb/build_module.py",
        "svdb/overlap_module.py",
        "svdb/readVCF.py",
        "svdb/merge_vcf_module_cython.py",
        "svdb/hist_module.py",
        "svdb/query_module.py",
        "svdb/export_module.py"])
else:
    ext_modules = []

setup(
    name = 'svdb',
    version = '1.0.2',
    ext_modules = ext_modules,
    packages = ['svdb'],
    install_requires = ['numpy', 'scikit-learn==0.15.2', 'scipy'],
    entry_points = {'console_scripts': ['svdb = svdb.__main__:main']},
)
