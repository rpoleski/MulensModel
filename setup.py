from setuptools import setup
from setuptools.command.install import install
import subprocess
import os
import sys


dir_ = os.path.join('lib', 'python' + sys.version[:3], 'site-packages', 
                    'MulensModel')

source_VBBL = os.path.join('source', 'VBBL')
source_AC = os.path.join('source', 'AdaptiveContouring')
source_MM = os.path.join('source', 'MulensModel')
source_MMmo = os.path.join(source_MM, 'mulensobjects')

wrapper_VBBL = os.path.join(source_VBBL, 'VBBinaryLensingLibrary_wrapper.so')
wrapper_AC = os.path.join(source_AC, 'AdaptiveContouring_wrapper.so')


class CustomInstall(install):
    """
    Custom install procedure that runs makefiles.
    """
    def run(self):
        print("Begin running makefiles...")
        kwargs = {"stderr": subprocess.STDOUT, 'shell': True}
        make_process = subprocess.Popen("make -C " + source_VBBL, **kwargs)
        if make_process.wait():
            print("VBBL MAKEFILE FAILED.")
        make_process = subprocess.Popen("make -C " + source_AC, **kwargs)
        if make_process.wait():
            print("AdaptiveContouring MAKEFILE FAILED.")
        print("Finish running makefiles...")
        install.run(self)

setup(
    name='MulensModel',
    version='1.2.10',
    url='git@github.com:rpoleski/MulensModel.git',
    cmdclass={'install': CustomInstall},
    author='Radek Poleski',
    author_email='poleski.1@osu.edu',
    description='packge for modeling gravitational microlensing events',
    packages=['MulensModel', 'MulensModel.mulensobjects'],
    package_dir = {
        'MulensModel': source_MM,
        'MulensModel.mulensobjects': source_MMmo},
    data_files=[
        (os.path.join(dir_, source_VBBL), [wrapper_VBBL]),
        (os.path.join(dir_, source_AC), [wrapper_AC])],
    install_requires=['numpy >= 1.11.1', 'matplotlib >= 1.5.1', 'scipy',
        'astropy']
)
