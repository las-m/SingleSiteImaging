# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 20:58:01 2024

@author: scott
"""

import add_path
from sisipy import sisisim
import numpy as np
import matplotlib.pyplot as plt

setup = sisisim.open_sim("setup50x50.pkl")

setup.changeEmitter('poisson', mu=200)

setup.sampleImages(500, plot=False, filling=0.2, save=True, savename='imgs50x50_200ph')