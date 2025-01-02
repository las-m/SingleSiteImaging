# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 08:45:47 2024

@author: scott
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
from sisi import quantum_emitters
from scipy.integrate import dblquad
import scipy.stats
import pickle

def open_sim(fname):
    with open(fname, "rb") as f:
        return pickle.load(f)

def gen_atom_plane(*lattice_params, **lattice_kwargs):
    """
    Checks an input as either an AtomPlane or a list of parameters to pass to AtomPlane __init__ function, then always returns an AtomPlane.

    Parameters
    ----------
    lattice_params : Lattice or list
        list of params
        
    lattice_kwargs: dict
        kwargs to pass to Lattice __init__

    Returns
    -------
    atom_plane : AtomPlane
        DESCRIPTION.

    """
    if type(lattice_params) == Lattice:
        return lattice_params
    elif lattice_params[0] is None:
        return Lattice(**lattice_kwargs)
    elif type(lattice_params) == list:
        return Lattice(*lattice_params, **lattice_kwargs)
    else:
        raise ValueError("lattice_params must be AtomPlane, list, or None")
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

class Lattice():
    def __init__(self, n_sites, spacing, *emitter_args, **emitter_kwargs):
        """
        Initialize a lattice structure with associated quantum emitters.

        Parameters
        ----------
        n_sites : int or array-like of int
            Number of lattice sites in each direction. If an integer is provided, 
            it will be used for all directions. If an array-like object is provided, 
            it should specify the number of lattice sites in each direction.
        
        spacing : float or array-like of float
            Lattice spacing along each direction (m). If a single float is provided, 
            it will be used for all directions. If an array-like object is provided, 
            it should specify the lattice spacing in each direction.

        pdist : str, array-like, or rv_discrete
            Probability distribution for the emitter photon counts. Can be a string
            for predefined distributions (e.g., 'poisson'), an array-like object
            for custom discrete probabilities, or an `rv_discrete` instance.

        *emitter_args : tuple
            Additional positional arguments for the `Emitter` class.

        **emitter_kwargs : dict
            Additional keyword arguments for the `Emitter` class.

        Raises
        ------
        ValueError
            If the dimensions of `n_sites` and `spacing` are incompatible.

        """
        n_sites_arr = np.atleast_1d(n_sites).astype(int).flatten()
        spacing_arr = np.atleast_1d(spacing).astype(float).flatten()
        self.emitter = quantum_emitters.Emitter(*emitter_args, **emitter_kwargs)

        if len(n_sites_arr) == len(spacing_arr) and len(n_sites_arr) != 1:
            self.n_sites = n_sites_arr
            self.spacing = spacing_arr
        elif len(n_sites_arr) == 1 and len(spacing_arr) == 1:
            self.n_sites = np.array([n_sites_arr,n_sites_arr]).flatten()
            self.spacing = np.array([spacing_arr,spacing_arr]).flatten()
        elif len(n_sites_arr) == 1 or len(spacing_arr) == 1:
            n_sites_arr, spacing_arr = np.meshgrid(n_sites_arr, spacing_arr)
            self.n_sites = n_sites_arr.flatten()
            self.spacing = spacing_arr.flatten()
        else:
            raise ValueError("Invalid dimensions")

        self.dim = len(self.n_sites)
        self.refreshCoordinates()


    def refreshCoordinates(self):
        x = [np.arange(0, self.n_sites[i])*self.spacing[i]-(self.n_sites[i]-1)*self.spacing[i]/2 for i in range(self.dim)]
        x_mesh = np.meshgrid(*x)
        self.coords = [x_mesh[i].flatten() for i in range(self.dim)]
        return
    
    def genDistribution(self, filling=1):
        return np.array([np.random.random(self.n_sites) < filling], dtype=int).flatten()
    
    def changeEmitter(self, *args, **kwargs):
        self.emitter = quantum_emitters.Emitter(*args, **kwargs)

class PSFSampler():
    def __init__(self, psf, psf_params):
        if callable(psf):
            self.psf = psf
        else:
            raise ValueError("psf must be callable")
        self.psf_params = psf_params
        
    def __call__(self, x, y):
        xy = x,y
        return self.psf(xy, *self.psf_params)
    
    # def integrate(self, x_min, x_max, y_min, y_max):
    #     integral, error = dblquad(self, x_min, x_max, lambda x: y_min, lambda x: y_max)
    #     return integral
    
    def integrate(self, x_min, x_max, y_min, y_max):
        if self((x_min+x_max)/2,(y_min+y_max)/2) > 1e-6:
            x, dx = np.linspace(x_min, x_max, 11, retstep=True)
            y, dy = np.linspace(y_min, y_max, 11, retstep=True)
            integral = np.sum(self(x,y))*dx*dy
        else:
            integral = 0
        return integral
    
class ImagingPlane():
    def __init__(self, n_pixels, pixel_size, magnification, psf, psf_params):
        if len([n_pixels]) == 1:
            self.n_pixels = np.array([n_pixels]*2, dtype=int).flatten()
        elif len([n_pixels]) == 2:
            self.n_pixels = np.array([n_pixels], dtype=int).flatten()
        else:
            raise ValueError('n_pixels must be int or array-like of length 2')
        
        self.pixel_size = float(pixel_size)
        self.magnification = float(magnification)
        self.psf = PSFSampler(psf, psf_params)
        self.psf_list = []
        
    def pixelate(self, x_atom, y_atom):
        """
        Returns a point-spread function for an atom centered at xy, by integrating
        the PSF over each camera pixel.

        Parameters
        ----------
        x : float
        
        y : float

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        
        psf_pixelated = np.zeros(self.n_pixels)
        
        for i in range(self.n_pixels[0]):
            for j in range(self.n_pixels[1]):
                x_min = (i - self.n_pixels[0]/2)*self.pixel_size/self.magnification - x_atom
                x_max = x_min + self.pixel_size/self.magnification
                y_min = (j - self.n_pixels[1]/2)*self.pixel_size/self.magnification - y_atom
                y_max = y_min + self.pixel_size/self.magnification
                
                psf_pixelated[i,j] = self.psf.integrate(x_min, x_max, y_min, y_max)
        
        psf_pixelated /= psf_pixelated.sum()
        
        return psf_pixelated
        
    def updatePSFList(self, x, y):
        """
        

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if len(x) == len(y):
            self.psf_list = [self.pixelate(x[i], y[i]) for i in range(len(x))]
        else:
            raise ValueError("Arrays must be the same length")
            
    def samplePSF(self, i, n_counts):
        """
        

        Parameters
        ----------
        i : int
            Index of PSF to sample.
            
        n_counts : int
            Number of times to sample.

        Returns
        -------
        None.

        """
        # Flatten the point mass for sampling
        indices = np.arange(self.n_pixels.prod())
        
        # Sample indices with the given weights
        sampled_indices = np.random.choice(indices, size=n_counts, p=self.psf_list[i].flatten())
        
        # Map sampled indices back to 2D grid
        image = np.zeros(self.n_pixels, dtype=int)
        for idx in sampled_indices:
            x, y = np.unravel_index(idx, self.n_pixels)
            image[x, y] += 1
        return image
    
    def sampleNoise(self, mu):
        pdist = scipy.stats.poisson(mu=mu)
        return pdist.rvs(size=self.n_pixels)
        
class Experiment(Lattice):
    def __init__(self, n_sites = None, spacing = None, **kwargs):
        super().__init__(n_sites, spacing, **kwargs)
        self.imaging_planes = []
        
    def addImagingPlane(self, imaging_plane=None, **kwargs):
        self.imaging_planes += [gen_imaging_plane(imaging_plane, **kwargs)]
        self.calcPSFs()
        
    def calcPSFs(self):
        for imaging_plane in self.imaging_planes:
            imaging_plane.updatePSFList(self.coords[0], self.coords[1])
        
    def sampleImages(self, n, filling=1, save=False, plot=False, snr_inf=False, noise=True):
        if save:
            print("Saving not implemented yet.")
        for i in range(n):
            atom_dist = self.genDistribution(filling).flatten()
            for imaging_plane in self.imaging_planes:
                img = np.zeros(imaging_plane.n_pixels)
                for k in range(self.n_sites.prod()):
                    if atom_dist[k]:
                        if snr_inf:
                            img += imaging_plane.psf_list[k]
                        else:
                            n_photons = self.emitter.sample_counts()
                            img += imaging_plane.samplePSF(k, n_photons)
                if noise and not snr_inf:
                    img += imaging_plane.sampleNoise(3)
                
                if plot:
                    plt.figure()
                    plt.imshow(img, cmap='RdBu_r', vmax=10)
                    plt.colorbar()
                    plt.show()
                    
    def save(self, fname):
        with open("%s.pkl" %fname, "wb") as file:
            pickle.dump(self, file)