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
n_sites = 20 # n x n sites

# Imaging parameters
magnification = 50
n_pixels = 200
pixel_size = 6.5e-6
psf_func = sisi.psf_functions.gaussian_2d
psf_params = [0, 0, 313e-9, 313e-9, 0]

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
    
    'psf'           : sisi.psf_functions.gaussian_2d,
    'psf_params'    : [1, 0, 0, 313e-9, 313e-9, 0] #amp, x_0, y_0, w_x, w_y, z0
    }

imaging_system = sisisim.Experiment(**atom_plane_params)
imaging_system.addImagingPlane(**imaging_plane_params)
imaging_system.sampleImages(1, save=False, plot=True, snr_inf=True)

#%%

# psf_disc = imaging_system.imaging_planes[0].psf_list
# psf_call = imaging_system.imaging_planes[0].psf

# x,y = np.meshgrid(np.linspace(-10e-6, 10e-6,100), np.linspace(-10e-6, 10e-6,100))
# plt.figure()
# plt.imshow(psf_call(x,y))
# plt.colorbar()
