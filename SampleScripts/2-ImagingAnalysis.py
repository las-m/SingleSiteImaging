# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 20:58:01 2024

@author: scott
"""

import add_path
from sisipy import sisisim, sisirecon
import numpy as np
import matplotlib.pyplot as plt

# fname = "imgs50x50_200ph.tiff"
fname = 'imgs50x50_1000ph.tiff'
setup = sisisim.open_sim("setup50x50.pkl")

imgs = sisirecon.ImageLoader(fname)

#%%

n_imgs = 1
n_avgs = 500
k = 0.172789251

ylim=None
# ylim=(50,3000)

[imgs.display(i, show_sites=True, cmap='RdBu_r') for i in range(n_imgs)]

pixel_sum_list = []
for i in range(0,500):
    pixel_sum_list += [n for n in imgs.getPixelSumList(i,3, prefilter=True) if n is not None]
# imgs.getHistogram(i,kernel_size=7,plot=False)

#%%
plt.figure()
plt.hist(pixel_sum_list, bins=99)

# [imgs.fourierFixedFreq(i, k, plot_mag=True, plot_phase=True, ylim=ylim) for i in range(n_imgs)]

# [imgs.fourierFixedFreq(i, k, plot=True, ylim=(50,3000)) for i in range(3)]

# for i in range(1):
#     f_x, f_y = sisirecon.get_frequency_array(imgs.imgs[i].shape)

#     theta=np.linspace(0,2*np.pi,1001)
#     x,y = np.cos(theta)*k,np.sin(theta)*k

#     plt.figure()
#     plt.imshow(np.log10(imgs.fourierTransform(i)[0]), extent=(f_x.min(),f_x.max(),f_y.min(),f_y.max()),cmap='RdBu_r', vmin=np.log10(50),vmax=np.log10(3000))
#     plt.plot(x,y, alpha=0.2, color='k', lw=4)
#     plt.colorbar()
#     plt.show()

#     plt.figure()
#     plt.imshow(imgs.fourierTransform(i)[1], extent=(f_x.min(),f_x.max(),f_y.min(),f_y.max()),cmap='RdBu_r')
#     # plt.plot(x,y, alpha=0.2, color='k', lw=4)
#     plt.colorbar()
#     plt.show()

# ft_mag = np.average([imgs.fourierFixedFreq(i, k)[1] for i in range(n_avgs)],axis=0)
# ft_phase = np.average([imgs.fourierFixedFreq(i, k)[2] for i in range(n_avgs)],axis=0)
# theta = np.linspace(0,2*np.pi,500, endpoint=False)

# plt.figure()
# plt.plot(theta, ft_mag)
# plt.show()

# plt.figure()
# plt.plot(theta, ft_phase)
# plt.show()

# plt.figure()
# plt.imshow(np.average([imgs.imgs[i] for i in range(500)], axis=0))
# plt.colorbar()

# plt.figure()
# plt.imshow(np.average([np.log10(imgs.fourierTransform(i)[0]) for i in range(n_avgs)], axis=0), vmax=2.5)
# plt.colorbar()

# plt.figure()
# plt.imshow(np.average([imgs.fourierTransform(i)[1] for i in range(n_avgs)], axis=0))
# plt.colorbar()