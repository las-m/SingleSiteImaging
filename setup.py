# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 21:10:20 2024

@author: scott
"""

from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Simulations and analysis of images from K39 quantum gas microscope'
LONG_DESCRIPTION = ''

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="singlesite", 
        version=VERSION,
        author="Scott Hubele",
        author_email="<scott.hubele@uni-hamburg.de>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[]
)
