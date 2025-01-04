# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 20:58:01 2024

@author: scott
"""

import add_path
from sisipy import sisisim, sisirecon
import numpy as np
import matplotlib.pyplot as plt

fname = "imgs50x50_1000ph.tiff"
setup = sisisim.open_sim("setup50x50.pkl")

imgs = sisirecon.ImageLoader(fname)

#%%

# [imgs.display(i, cmap='RdBu_r') for i in range(3)]

for a in np.linspace(0.5,1.4,9):
    k = 0.172789251*a
    
    # [imgs.fourierFixedFreq(i, k, plot=True, ylim=(50,3000)) for i in range(3)]
    
    for i in range(1):
        f_x, f_y = sisirecon.get_frequency_array(imgs.imgs[i].shape)
        
        theta=np.linspace(0,2*np.pi,1001)
        x,y = np.cos(theta)*k,np.sin(theta)*k
        
        
        plt.figure()
        plt.imshow(np.log10(imgs.fourierTransform(i)), extent=(f_x.min(),f_x.max(),f_y.min(),f_y.max()),cmap='RdBu_r', vmin=np.log10(50),vmax=np.log10(3000))
        plt.plot(x,y, alpha=0.2, color='k', lw=4)
        plt.colorbar()
    
    ft = np.average([imgs.fourierFixedFreq(i, k)[1] for i in range(500)],axis=0)
    theta = np.linspace(0,2*np.pi,500, endpoint=False)
    
    plt.figure()
    plt.plot(theta, ft)
    plt.show()

# plt.figure()
# plt.imshow(np.average([imgs.imgs[i] for i in range(500)], axis=0))
# plt.colorbar()

# plt.figure()
# plt.imshow(np.average([imgs.fourierTransform(i) for i in range(500)], axis=0))
# plt.colorbar()