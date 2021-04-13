from distutils.core import setup
import setuptools
setup(name='PYT-SBI',
    version='1.0',
    description='PYT-SBI project',
    author='Mireia Codina Tobías, Ángel Monsalve Fernández, María Torralvo Márquez',
    url='https://github.com/angelmf97/PYT-project',
    packages=setuptools.find_packages(where='.')    # Change this in case we finally put the packages inside a directory
    )