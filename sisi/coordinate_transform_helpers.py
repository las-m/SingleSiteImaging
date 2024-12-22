# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 08:45:47 2024

@author: scott
"""

import numpy as np
import matplotlib.pyplot

class PSF():
    def __init__(self, xy, PSF_discreet):
        self.xy = xy
        self.PSF_discreet = PSF_discreet

class AtomPlane():
    def __init__(self, fov, spacing, origin=(0,0)):
        self.fov = float(fov)
        self.spacing = float(spacing)
        if len(origin)==2:
            self.origin=(float(origin[0]),float(origin[1]))
        else:
            raise ValueError("Origin must be (x,y) format")
            
class ImagingPlane():
    def __init__(self, magnification, pixel_size, origin=(0,0), centered="origin"):
        self.magnification=magnification
        if centered == "corners":    
            self.centered = "corners"
        elif centered == "origin":
            self.centered = "origin"
        else:
            raise ValueError("""Must be either "corners"-centered or "origin"-centered""")
        if len(origin)==2:
            self.origin=(float(origin[0]),float(origin[1]))
        else:
            raise ValueError("Origin must be (x,y) format")
    
class Experiment():
    def __init__(self, AtomPlane, *ImagingPlanes):
        self.AtomPlane = AtomPlane
        self.ImagingPlanes = [*ImagingPlanes]
        
    def addImagingPlane(self, *ImagingPlanes):
        self.ImagingPlanes+=[*ImagingPlanes]
        