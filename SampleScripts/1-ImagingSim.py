# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 20:58:01 2024

@author: scott
"""

import add_path
from sisi import sisisim
import sisi.psf_functions
import numpy as np
import matplotlib.pyplot as plt

# Atom plane parameters
n_sites = 50 # n x n sites

# Imaging parameters
magnification = 50
n_pixels = 200
pixel_size = 6.5e-6
psf_func = sisi.psf_functions.gaussian_2d
psf_params = [1, 313e-9, 313e-9] #amp, w_x, w_y

atom_plane_params = {
    'n_sites'    : n_sites,
    'spacing'    : 1064e-9/np.sqrt(2),
    'pdist'      : 'poisson',
    'mu'         : 10000
    }

imaging_plane_params = {
    'magnification' : magnification,
    'pixel_size'    : pixel_size,
    'n_pixels'      : n_pixels,
    
    'psf'           : psf_func,
    'psf_params'    : psf_params,
    'extent'        : 1e-5
    }

imaging_system = sisisim.Experiment(**atom_plane_params)
imaging_system.addImagingPlane(**imaging_plane_params)
imaging_system.save("setup50x50")
imaging_system.sampleImages(1, save=False, plot=True, snr_inf=True)