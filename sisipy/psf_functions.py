# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 22:18:35 2024

List of possible PSFs. Format should be:
    func(xy, *args)
    
    where xy = x,y are the coordinates in the atom-plane in meters.

PSFs should have their center at (0,0)

@author: scott
"""

import numpy as np

def gaussian_2d(xy, a, waist_x, waist_y):
    """
    Compute the value of a 2D Gaussian function at a given point (x, y).

    Parameters:
        
        xy : xy 
            A tuple (x, y) representing the point where the Gaussian is evaluated.
        a : float
            Amplitude (peak value) of Gaussian.
        waist_x : float 
            The standard deviation of the Gaussian along the x-axis.
        waist_y : float 
            The standard deviation of the Gaussian along the y-axis.

    Returns:
        float: The computed value of the 2D Gaussian at the given point (x, y).
    """
    x, y = xy
    exponent = -((x ** 2) / (2 * waist_x ** 2) + (y ** 2) / (2 * waist_y ** 2))
    return a * np.exp(exponent)