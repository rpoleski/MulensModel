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

# from distutils.core import Extension
# ext_module = Extension("some_NAME_HERE", 
#   sources=["..."], # relative to the distribution root (where the setup script lives)
#   include_dirs=["...],
#   library_dirs = [os.getcwd(),],  # path to .a or .so file(s)
#   extra_compile_args=["-O2", "-lm"],
#   extra_link_args=["", ""], # this will be for last command only, most probably
#   depends=["", ""], # list of files that the extension depends on
#   optional=False/True, # specifies that a build failure in the extension should not abort the build process, but simply skip the extension.
#   language='c++11']
#later:
#    ext_modules=[ext_module]

# import subprocess
# make_process = subprocess.Popen("make clean all", stderr=subprocess.STDOUT)
# if make_process.wait() != 0:
#     DO_SOMETHING

#from distutils.command.install import install
# OR:
#from setuptools.command.install import instal
#
#class CustomInstall(install):
#def run(self):
#    install.run(self)
#    # custom stuff here
#    do_my_stuff()
#
#setup(..., cmdclass={'install': CustomInstall})

# Generate documentation dictionary and save it in "lib/"
# import get_docstring
# docstring = get_docstring.get_docstring()
# f = open("lib/docstring_pickle.pkl", "wb")
# pickle.dump(docstring, f)
# f.close()
#and later setup options were:
#package_dir = {'pyslalib': 'lib'},
#package_data = {'pyslalib': ['docstring_pickle.pkl']},
