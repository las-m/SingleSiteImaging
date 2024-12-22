# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 04:34:14 2024

@author: scott
"""

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian_2d(xy, A, mu_x, mu_y, sigma_x, sigma_y, z0):
    x, y = xy
    return z0 + A * np.exp(-((x - mu_x)**2 / (2 * sigma_x**2) + (y - mu_y)**2 / (2 * sigma_y**2)))

pixel_size = 6.5e-6/50 # Equivalent pixel size at atom plane in um

x = np.linspace(-99.5*pixel_size, 99.5*pixel_size, 200)
y = np.linspace(-99.5*pixel_size, 99.5*pixel_size, 200)
Xmesh,Ymesh = np.meshgrid(x,y)

shots = [np.asarray(Image.open("SingleAtom_Shot%i.png" %i)) for i in range(1,6)]
shot_avg = np.average(shots, axis=0)

shot_avg = (shot_avg - shot_avg.min())/(shot_avg.max()-shot_avg.min())
argmax = [shot_avg.argmax()//200, shot_avg.argmax()%200]

avg_tightcrop = shot_avg[(argmax[0]-12):(argmax[0]+12),(argmax[1]-12):(argmax[1]+12)]
Xmesh_tightcrop = Xmesh[(argmax[0]-12):(argmax[0]+12),(argmax[1]-12):(argmax[1]+12)]
Ymesh_tightcrop = Ymesh[(argmax[0]-12):(argmax[0]+12),(argmax[1]-12):(argmax[1]+12)]

#%%

xfit = np.linspace(Xmesh_tightcrop.min(), Xmesh_tightcrop.max(), 1001)
yfit = np.linspace(Ymesh_tightcrop.min(), Ymesh_tightcrop.max(), 1001)

XfitMesh, YfitMesh = np.meshgrid(xfit,yfit)

initial_guess = (1, x[argmax[1]], y[argmax[0]], 1e-6, 1e-6, 0)
params,_ = curve_fit(gaussian_2d, (Xmesh_tightcrop.ravel(), Ymesh_tightcrop.ravel()), avg_tightcrop.ravel(), p0=initial_guess)

# Plot the fitted Gaussian
psf_fit = gaussian_2d((Xmesh, Ymesh), *params)
psf_fit_crop = gaussian_2d((XfitMesh, YfitMesh), *params)

#%%

extent = np.array([-(xfit.max()-xfit.min())/2, (xfit.max()-xfit.min())/2, -(yfit.max()-yfit.min())/2, (yfit.max()-yfit.min())/2])/1e-6

# Plot the original data and fitted data side by side
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

im00 = axes[0][0].imshow(shot_avg, vmin=0, vmax=1)
axes[0][0].set_title("Average of 5 Shots")
axes[0][0].set_xlabel('x (px)')
axes[0][0].set_ylabel('y (px)')

im01 = axes[0][1].imshow(psf_fit, extent=extent, vmin=0, vmax=1)
axes[0][1].set_title("Gaussian fit - uncropped")
axes[0][1].set_xlabel('x (px)')
axes[0][1].set_ylabel('y (px)')

# Original data plot
im10 = axes[1][0].imshow(avg_tightcrop, extent=extent, vmin=0, vmax=1)
axes[1][0].set_title("Average shots - cropped")
axes[1][0].set_xlabel('x (μm)')
axes[1][0].set_ylabel('y (μm)')

# Fitted data plot
im11 = axes[1][1].imshow(psf_fit_crop, extent=extent, vmin=0, vmax=1)
axes[1][1].set_title("Gaussian fit - cropped")
axes[1][1].set_xlabel('x (μm)')
axes[1][1].set_ylabel('y (μm)')

# Add a shared colorbar
fig.colorbar(im11, ax=axes, orientation='vertical', fraction=0.06, pad=0.1)

plt.savefig('GaussianFit.png')

#%%

# plt.figure()
# plt.imshow(shot_avg, vmin=0, vmax=1)
# plt.colorbar()
# plt.plot(argmax[1],argmax[0], 'ko')

# plt.figure()
# plt.imshow(avg_tightcrop)

# fig = plt.figure()
# plt.subplot(211)
# plt.plot(x, shot_avg[argmax[0],:])
# plt.subplot(212)
# plt.plot(y, shot_avg[:,argmax[1]])

