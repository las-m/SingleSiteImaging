# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 04:13:12 2025

@author: scott
"""

import imageio
import numpy as np
import matplotlib.pyplot as plt
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
    def __init__(self, fname):
        self.imgs = imageio.mimread(fname)
        
    def display(self, i, **kwargs):
        plt.figure()
        plt.imshow(self.imgs[i], **kwargs)
        plt.show()
        
    def fourierTransform(self, i, plot=False):
        fft_image = np.fft.fft2(self.imgs[i])
        fft_shifted = np.fft.fftshift(fft_image)
        magnitude_spectrum = np.abs(fft_shifted)
        return magnitude_spectrum
    
    def fourierFixedFreq(self, i, freq, n_pts = 500, plot=False, ylim=None, **kwargs):
        mag_spectrum = self.fourierTransform(i)
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
        
        if plot:
            fig, ax = plt.subplots()
            ax.semilogy(theta*180/np.pi, mag_fixed_freq, **kwargs)
            if ylim is not None:
                ax.set_ylim(ylim)
            ax.set_xlabel(r'Angle ($\deg$)')
            ax.set_ylabel('Fourier Amplitude')
            plt.show()
            
        return theta, mag_fixed_freq
    
    def getLatticePosition(self, i, method='ft', freq=0.17278925104182):
        theta, ft_amp = self.fourierFixedFreq(i, freq)