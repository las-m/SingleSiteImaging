# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 20:58:01 2024

@author: scott
"""

import add_path
from sisipy import sisisim
import sisipy.psf_functions
import numpy as np
import matplotlib.pyplot as plt

# Atom plane parameters
n_sites = 3 # n x n sites

# Imaging parameters
magnification = 50
n_pixels = 30
pixel_size = 6.5e-6
psf_func = sisipy.psf_functions.gaussian_2d
psf_params = [1, 313e-9, 313e-9] #amp, w_x, w_y

atom_plane_params = {
    'n_sites'    : n_sites,
    'spacing'    : 1064e-9/np.sqrt(2),
    'pdist'      : 'poisson',
    'mu'         : 100,
    'phase'      : (2*np.pi/3,2*np.pi/3)
    }

imaging_plane_params = {
    'magnification' : magnification,
    'pixel_size'    : pixel_size,
    'n_pixels'      : n_pixels,
    'angle'         : 10,
    
    'psf'           : psf_func,
    'psf_params'    : psf_params,
    'extent'        : 1e-5
    }

imaging_system = sisisim.Experiment(**atom_plane_params)
imaging_system.addImagingPlane(**imaging_plane_params)
# imaging_system.save("setup50x50")
imaging_system.sampleImages(1, save=False, filling=1, plot=True, snr_inf=True)
# imaging_system.sampleImages(500, plot=False, filling=0.2, save=True, savename='imgs50x50_100ph_120deg_phase')