# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 04:13:12 2025

@author: scott
"""

import imageio
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib.ticker import FuncFormatter

def get_frequency_array(image_shape, axes='cartesian'):
    # Get the dimensions of the image
    rows, cols = image_shape
    
    # Compute the frequency coordinates
    freq_y = np.fft.fftfreq(rows)  # Frequencies for rows
    freq_x = np.fft.fftfreq(cols)  # Frequencies for columns
    
    # Shift the zero frequency to the center
    freq_y = np.fft.fftshift(freq_y)
    freq_x = np.fft.fftshift(freq_x)
    
    if axes=='cartesian':
        return freq_x, freq_y
    elif axes=='radial':
        # Create 2D grids of frequencies
        f_x, f_y = np.meshgrid(freq_x, freq_y)
        
        # Compute the frequency magnitude
        f_r = np.sqrt(f_x**2 + f_y**2)
        angle = np.arctan(f_y/f_x)*180/np.pi
        return f_r, angle
    else:
        raise ValueError("axes must be either 'radial' or 'cartesian'")

class ImageLoader():
    def __init__(self, fname, setup=None):
        self.imgs = imageio.mimread(fname)
        
        self.pixel_size = 6.5e-6
        self.shape = np.array(self.imgs[0].shape)
        
    def display(self, i, prefilter=False, show_sites=False, **kwargs):
        if prefilter:
            img = self.gaussianFilter(i)
        else:
            img = self.imgs[i]
        plt.figure()
        plt.imshow(img, **kwargs)
        plt.xlim([0,self.shape[1]])
        plt.ylim([self.shape[0],0])
        if show_sites:
            origin = self.getLatticeOrigin(i)
            x, y = self.getLatticeCoords(i, origin, self.getLatticeAngle(i))
            plt.plot(x,y, 'ko', markersize=2)
    
    def gaussianFilter(self, i, sigma=0.5, **kwargs):
        return ndimage.gaussian_filter(self.imgs[i], sigma)
        
    def fourierTransform(self, i, plot=False, prefilter=False):
        if prefilter:
            fft_image = np.fft.fft2(self.gaussianFilter(i))
        else:
            fft_image = np.fft.fft2(self.imgs[i])
        fft_shifted = np.fft.fftshift(fft_image)
        magnitude_spectrum = np.abs(fft_shifted)
        phase_spectrum = np.angle(fft_shifted)
        return magnitude_spectrum, phase_spectrum
    
    def fourierFixedFreq(self, i, freq, n_pts = 500, plot_mag=False, plot_phase=False, ylim=None, prefilter=False,**kwargs):
        mag_spectrum, phase_spectrum = self.fourierTransform(i, prefilter=prefilter)
        theta = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
        
        freq_x = np.cos(theta)*freq
        freq_y = np.sin(theta)*freq
        
        n_rows, n_cols = self.imgs[i].shape
        
        freq_y_scale = np.fft.fftfreq(n_rows)  # Frequencies for rows
        freq_x_scale = np.fft.fftfreq(n_cols)  # Frequencies for columns
        
        freq_y_scale = np.fft.fftshift(freq_y_scale)
        freq_x_scale = np.fft.fftshift(freq_x_scale)
        
        freq_x_mat, freq_x_scale_mat = np.meshgrid(freq_x, freq_x_scale)
        freq_y_mat, freq_y_scale_mat = np.meshgrid(freq_y, freq_y_scale)
        
        closest_freq_idx = np.argmin(np.abs(freq_x_mat-freq_x_scale_mat), axis=0)
        closest_freq_idy = np.argmin(np.abs(freq_y_mat-freq_y_scale_mat), axis=0)
        
        mag_fixed_freq = np.array([mag_spectrum[closest_freq_idy[i], closest_freq_idx[i]] for i in range(n_pts)])
        phase_fixed_freq = np.array([phase_spectrum[closest_freq_idy[i], closest_freq_idx[i]] for i in range(n_pts)])
        
        if plot_mag:
            fig, ax = plt.subplots()
            ax.semilogy(theta*180/np.pi, mag_fixed_freq, **kwargs)
            if ylim is not None:
                ax.set_ylim(ylim)
            ax.set_xlabel(r'Angle ($\deg$)')
            ax.set_ylabel('Fourier Amplitude')
        
        if plot_phase:
            fig, ax = plt.subplots()
            ax.plot(theta*180/np.pi, phase_fixed_freq, **kwargs)
            if ylim is not None:
                ax.set_ylim(ylim)
            ax.set_xlabel(r'Angle ($\deg$)')
            ax.set_ylabel('Fourier Phase') 
        
        return theta, mag_fixed_freq, phase_fixed_freq
    
    def getLatticeOrigin(self, i, method='ft', k=None):
        if k is None:
            k=self.getLatticeFreq(i)
        theta, ft_amp, ft_phase = self.fourierFixedFreq(i, k)
        
        lat_phase = np.array([2.7,2.52]) # x,y / axis1,axis2
        
        # origin = self.pixel_size*((np.array(self.shape)/2)%k + lat_phase/(2*np.pi*k))
        origin = lat_phase/(2*np.pi*k)
        # origin = np.array([100,100])
        
        return origin
    
    def getLatticeAngle(self, i, k=None):
        if k is None:
            k=self.getLatticeFreq(i)
        return 10*np.pi/180
    
    def getLatticeFreq(self, i):
        return 0.17278925104182
    
    def getLatticeCoords(self, i, origin, theta, k=None):
        if k is None:
            k=self.getLatticeFreq(i)
        x_0, y_0 = origin
        n_x, n_y = 4*self.shape*k
        
        x = np.arange(-n_x//2, n_x//2)/k
        y = np.arange(-n_y//2, n_y//2)/k
        
        x,y = np.meshgrid(x,y)
        
        x_rot = np.cos(theta)*x + np.sin(theta)*y
        y_rot = -np.sin(theta)*x + np.cos(theta)*y
        
        x_rot += x_0
        y_rot += y_0
        
        in_bounds = np.all([(x_rot>0), (x_rot<self.shape[0]), (y_rot>0), (y_rot<self.shape[1])],axis=0)
        
        # plt.figure()
        # plt.plot(x.flatten(),y.flatten(),'ko',markersize=3)
        # plt.plot(x_rot.flatten(),y_rot.flatten(),'ro',markersize=3)
        
        return x_rot[in_bounds], y_rot[in_bounds]
        
    def sumPixels(self, i, xy, kernel_size=None, prefilter=False):

        if kernel_size is None:
            kernel_size = int(0.5//(self.getLatticeFreq(i)))
            kernel_size = 2*kernel_size + 1 #Ensures kernel size is odd
        elif int(kernel_size)%2 != 1:
            raise Warning("kernel_size must be an odd integer. Adjusting to odd kernel_size...")
            kernel_size = 2*(kernel_size//2)+1
    
        if prefilter:
            img = self.gaussianFilter(i)
        else:
            img = self.imgs[i]
    
        idx, idy = np.floor(xy).astype(int)
        if self.checkKernel(idx,idy,kernel_size):
            slice_x = slice(int(idx-kernel_size//2),int(idx+1+kernel_size//2))
            slice_y = slice(int(idx-kernel_size//2),int(idx+1+kernel_size//2))
            return np.sum(img[slice_x,slice_y])
        else:
            return None
        
    def checkKernel(self, idx, idy, kernel_size):
        x_is_in = ((idx-(kernel_size//2))>=0) and ((idx+1+kernel_size//2)<self.shape[0])
        y_is_in = ((idy-(kernel_size//2))>=0) and ((idy+1+kernel_size//2)<self.shape[1])
        return (x_is_in and y_is_in)
    
    def getPixelSumList(self, i, kernel_size=None, prefilter=True):
        site_pos_x, site_pos_y = self.getLatticeCoords(i, self.getLatticeOrigin(i), self.getLatticeAngle(i))
        
        pixel_sum_list = []
        
        for idx in range(len(site_pos_x)):
            x,y = site_pos_x[idx], site_pos_y[idx]
            pixel_sum_list.append(self.sumPixels(i, (x,y), kernel_size,prefilter=prefilter))
            
        return pixel_sum_list
    
    def getHistogram(self, i, kernel_size=None, plot=False, prefilter=True):
        pixel_sum_list = [n for n in self.getPixelSumList(i,kernel_siz,prefilter=prefilter) if n is not None]
        if plot:
            plt.figure()
            plt.hist(pixel_sum_list,bins=24)
        return np.histogram(pixel_sum_list)