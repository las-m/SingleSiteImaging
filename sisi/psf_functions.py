# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 22:18:35 2024

@author: scott
"""

import numpy as np

def gaussian_2d(xy, a, center_x, center_y, waist_x, waist_y, z0):
    """
    Compute the value of a 2D Gaussian function at a given point (x, y).

    Parameters:
        
        xy : xy 
            A tuple (x, y) representing the point where the Gaussian is evaluated.
        a : float
            Amplitude (peak value) of Gaussian.
        center_x : float 
            The x-coordinate of the Gaussian's center.
        center_y : float 
            The y-coordinate of the Gaussian's center.
        waist_x : float 
            The standard deviation of the Gaussian along the x-axis.
        waist_y : float 
            The standard deviation of the Gaussian along the y-axis.
        z0 : float
            Offset

    Returns:
        float: The computed value of the 2D Gaussian at the given point (x, y).
    """
    x, y = xy
    exponent = -(((x - center_x) ** 2) / (2 * waist_x ** 2) + ((y - center_y) ** 2) / (2 * waist_y ** 2))
    return a * np.exp(exponent) + z0