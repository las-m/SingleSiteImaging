# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 20:58:01 2024

@author: scott
"""

import add_path
from sisi import sisisim
import numpy as np
import matplotlib.pyplot as plt

setup = sisisim.open_sim("setup1.pkl")

setup.changeEmitter('poisson', mu=200)

setup.sampleImages(1, plot=True, filling=0.1)