import os
import sys
import glob
import warnings
from setuptools import setup, Extension, find_packages


file_required = "requirements.txt"

source_VBBL_abs = os.path.abspath(os.path.join('source', 'VBBL'))
source_AC_abs = os.path.abspath(os.path.join('source', 'AdaptiveContouring'))
source_MM = os.path.join('source', 'MulensModel')
source_MMmo = os.path.join(source_MM, 'mulensobjects')

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
dir_files["."] = [file_required, "README.md", "LICENSE"]
files_ = ["Horizons_manual.md", "examples_list.md", "papers_to_cite.md",
          "magnification_methods.pdf", "parameter_names.pdf"]
dir_files["documents"] = [os.path.join("documents", f) for f in files_]
dir_files["examples"] = glob.glob(os.path.join("examples", "*.py"))
dir_ex_16 = os.path.join("examples", "example_16")
dir_ex_16_files = ["ulens_model_fit.py", "ulens_model_plot.py"]
dir_files[dir_ex_16] = [os.path.join(dir_ex_16, f) for f in dir_ex_16_files]
data_files = []
for (key, value) in dir_files.items():
    data_files += [(os.path.join('MulensModel', key), value)]

# Read ther version.
version = "unknown"
with open(os.path.join('source', 'MulensModel', 'version.py')) as in_put:
    for line_ in in_put.readlines():
        if line_.startswith('__version__'):
            version = line_.split()[2][1:-1]

# Prepare extensions and required modules.
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
#        'git clone': 'git@github.com:rpoleski/MulensModel.git',
        'documentation': 'https://github.com/rpoleski/MulensModel'},
    ext_modules=[ext_AC, ext_VBBL],
    author='Radek Poleski & Jennifer Yee',
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
