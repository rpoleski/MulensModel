from setuptools import setup, find_packages

setup(
    name='MulensModel',
    version='1.2.10',
    url='git@github.com:rpoleski/MulensModel.git',
    author='Radek Poleski',
    author_email='poleski.1@osu.edu',
    description='packge for modeling gravitational microlensing events',
    packages=['source/MulensModel'],
    install_requires=['numpy >= 1.11.1', 'matplotlib >= 1.5.1'],
)

