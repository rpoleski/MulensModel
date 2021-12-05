import os
import sys
import glob
import warnings
from setuptools import setup, Extension, find_packages


file_required = "requirements.txt"

source_VBBL = os.path.join('source', 'VBBL')
source_AC = os.path.join('source', 'AdaptiveContouring')
source_MM = os.path.join('source', 'MulensModel')
source_MMmo = os.path.join(source_MM, 'mulensobjects')

source_VBBL_abs = os.path.abspath(source_VBBL)
source_AC_abs = os.path.abspath(source_AC)

# Read all files from data/ in format adequate for data_files option of setup.
files = glob.glob(os.path.join("data", "**", "*"), recursive=True)
files = [f for f in files if os.path.isfile(f)]
dir_files = dict()
for file_ in files:
    dir_ = os.path.dirname(file_)
    if dir_ in dir_files:
        dir_files[dir_] += [file_]
    else:
        dir_files[dir_] = [file_]
data_files = []
dir_files["."] = [file_required]
for (key, value) in dir_files.items():
    data_files += [(os.path.join('MulensModel', key), value)]

version = "unknown"
with open(os.path.join('source', 'MulensModel', 'version.py')) as in_put:
    for line_ in in_put.readlines():
        if line_.startswith('__version__'):
            version = line_.split()[2][1:-1]

ext_AC = Extension('MulensModel.AdaptiveContouring', libraries=["m"],
                   sources=glob.glob(os.path.join(source_AC_abs, "*.c")))
ext_VBBL = Extension('MulensModel.VBBL', libraries=["m"],
                     sources=glob.glob(os.path.join(source_VBBL_abs, "*.cpp")))

with open(file_required) as file_:
    required = file_.read().splitlines()

setup(
    name='MulensModel',
    version=version,
    url='https://github.com/rpoleski/MulensModel',
    project_urls={
        'git clone': 'git@github.com:rpoleski/MulensModel.git',
        'documentation': 'https://github.com/rpoleski/MulensModel'},
    ext_modules=[ext_AC, ext_VBBL],
    author='Radek Poleski',
    author_email='radek.poleski@gmail.com',
    description='package for modeling gravitational microlensing events',
    packages=['MulensModel', 'MulensModel.mulensobjects'],
    # This also should work: packages=find_packages(where="source"),
    package_dir={
        'MulensModel': source_MM,
        'MulensModel.mulensobjects': source_MMmo},
    data_files=data_files,
    install_requires=required,
    python_requires=">=3.6",
)
