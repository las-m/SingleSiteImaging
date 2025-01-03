# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 21:10:20 2024

@author: Scott Hubele
"""

from setuptools import setup, find_packages

VERSION = '0.0.3' 
DESCRIPTION = 'Simulations and analysis of images from K39 quantum gas microscope'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="sisi", 
        version=VERSION,
        author="Scott Hubele",
        url='git@github.com/las-m/SingleSiteImaging.git',
        author_email="<scott.hubele@uni-hamburg.de>",
        description=DESCRIPTION,
        packages=find_packages(),
        install_requires=[]
)
