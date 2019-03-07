import os
import sys
import subprocess
import glob
import warnings
from setuptools import setup, Extension


dir_ = os.path.join('lib', 'python' + sys.version[:3], 'site-packages', 
                    'MulensModel')

source_VBBL = os.path.join('source', 'VBBL')
source_AC = os.path.join('source', 'AdaptiveContouring')
source_MM = os.path.join('source', 'MulensModel')
source_MMmo = os.path.join(source_MM, 'mulensobjects')

data_files = [ (os.path.join(dir_, 'data'), [os.path.join('data', 'interpolation_table_b0b1_v1.dat')]) ]

version = "unknown"
with open(os.path.join('source', 'MulensModel', 'version.py')) as in_put:
    for line_ in in_put.readlines():
        if line_.startswith('__version__'):
            version = line_.split()[2][1:-1]

ext_AC = Extension('MulensModel.AdaptiveContouring', glob.glob(source_AC+"/*.c"))
ext_VBBL = Extension('MulensModel.VBBL', glob.glob(source_VBBL+"/*.cpp"))

setup(
    name='MulensModel',
    version=version,
    url='git@github.com:rpoleski/MulensModel.git',
    ext_modules=[ext_AC, ext_VBBL],
    author='Radek Poleski',
    author_email='poleski.1@osu.edu',
    description='packge for modeling gravitational microlensing events',
    packages=['MulensModel', 'MulensModel.mulensobjects'],
    package_dir={
        'MulensModel': source_MM,
        'MulensModel.mulensobjects': source_MMmo},
    data_files=data_files,
    install_requires=['numpy', 'matplotlib', 'scipy', 'astropy>=1.2.0']
)

