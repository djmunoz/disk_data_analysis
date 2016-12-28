from __future__ import print_function
from distutils.core import setup
import os
import sys


setup(
    name='disk_data_analysis',
    version="0.0.1",
    author='Diego J. Munoz',
    author_email = 'diego.munoz.anguita@gmail.com',
    url='https://github.com/',
    packages=['disk_data_analysis','disk_data_analysis.circumbinary'],
    description='Analysis and Plotting Tools for HDF5 simulation data of meshless/moving-mesh hydrodynamical disks',
    install_requires = ['numpy','scipy']
)
