# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 08:45:47 2024

@author: scott
"""

import numpy as np
import matplotlib.pyplot

def gen_atom_plane(*atom_plane_params, **atom_plane_kwargs):
    """
    Checks an input as either an AtomPlane or a list of parameters to pass to AtomPlane __init__ function, then always returns an AtomPlane.

    Parameters
    ----------
    atom_plane_params : AtomPlane or list
        list of params
        
    atom_plane_kwargs: dict
        kwargs to pass to AtomPlane __init__

    Returns
    -------
    atom_plane : AtomPlane
        DESCRIPTION.

    """
    if type(atom_plane_params) == AtomPlane:
        return atom_plane_params
    elif atom_plane_params[0] is None:
        return AtomPlane(**atom_plane_kwargs)
    elif type(atom_plane_params) == list:
        return AtomPlane(*atom_plane_params, **atom_plane_kwargs)
    else:
        raise ValueError("atom_plane_params must be AtomPlane, list, or None")
    return        

def gen_imaging_plane(*imaging_plane_params, **imaging_plane_kwargs):
    """
    Checks an input as either an ImagingPlane or a list of parameters to pass to ImagingPlane __init__ function, then always returns an ImagingPlane.

    Parameters
    ----------
    imaging_plane_params : ImagingPlane or list
        list of params
        
    imaging_plane_kwargs: dict
        kwargs to pass to ImagingPlane __init__

    Returns
    -------
    imaging_plane : ImagingPlane
        DESCRIPTION.

    """
    if type(imaging_plane_params) == ImagingPlane:
        return imaging_plane_params
    elif imaging_plane_params[0] is None:
        return ImagingPlane(**imaging_plane_kwargs)
    elif type(imaging_plane_params) == list:
        return ImagingPlane(*imaging_plane_params, **imaging_plane_kwargs)
    else:
        raise ValueError("imaging_plane_params must be ImagingPlane, list, or None")

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
    def __init__(self, AtomPlane=None, *ImagingPlanes):
        self.AtomPlane = AtomPlane
        self.ImagingPlanes = [*ImagingPlanes]
        
    def addImagingPlane(self, imaging_plane=None, **kwargs):
        self.ImagingPlanes+=[gen_imaging_plane(imaging_plane, **kwargs)]
        
    def setAtomPlaneParams(self, atom_plane=None, **kwargs):
        self.AtomPlane = gen_atom_plane(atom_plane, **kwargs)
        