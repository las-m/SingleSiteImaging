# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 20:58:01 2024

@author: scott
"""

import add_path
import sisi.coordinate_transform_helpers as sisisim
import sisi.psf_functions
import numpy as np

# Atom plane parameters
fov = 30e-6

# Imaging parameters
magnification = 50
pixel_size = 6.5e-6
psf_func = sisi.psf_functions.gaussian_2d
psf_params = [0, 0, 313e-9, 313e-9, 0]

atom_plane_params = {
    'fov'        : fov,
    'spacing'    : 1064e-9/np.sqrt(2)
    }

imaging_plane_params = {
    'magnification' : magnification,
    'pixel_size'    : pixel_size
    # 'psf_func'      : psf_func,
    # 'psf_params'    : psf_params
    }

imaging_system = sisisim.Experiment()
imaging_system.setAtomPlaneParams(**atom_plane_params)
imaging_system.addImagingPlane(**imaging_plane_params)