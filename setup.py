from setuptools import setup

setup(
    name='MulensModel',
    version='1.2.10',
    url='git@github.com:rpoleski/MulensModel.git',
    author='Radek Poleski',
    author_email='poleski.1@osu.edu',
    description='packge for modeling gravitational microlensing events',
    package_dir = {
        'MulensModel': 'source/MulensModel',
        'MulensModel.mulensobjects': 'source/MulensModel/mulensobjects'},
    packages=['MulensModel', 'MulensModel.mulensobjects'],
    data_files=[
        ('source/VBBL', ['source/VBBL/VBBinaryLensingLibrary_wrapper.so']),
        ('source/AdaptiveContouring', ['source/AdaptiveContouring/AdaptiveContouring_wrapper.so'])],
    install_requires=['numpy >= 1.11.1', 'matplotlib >= 1.5.1', 'math', 
        'os', 'sys', 'warnings', 'ctypes', 'scipy', 'astropy'],
)

# TO DO:
# - run makefile
# - 2 vs 3 levels in binarylens.py - probably we have to remove first "source/" in data_files records above
# - b0b1 file
# - install_requires
# - other setup options

# import subprocess
# make_process = subprocess.Popen("make clean all", stderr=subprocess.STDOUT)
# if make_process.wait() != 0:
#     DO_SOMETHING
