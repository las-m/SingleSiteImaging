# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 07:32:13 2024

@author: scott
"""

import numpy as np
import scipy.stats

class Emitter():
    def __init__(self, pdist='poisson', *args, **kwargs):
        if pdist == 'poisson' or pdist == 'Poisson':
            self.pdist = scipy.stats.poisson(*args, **kwargs)
        elif isinstance(pdist, (list, tuple, np.ndarray)):
            self.pdist = scipy.stats.rv_discrete(values=(np.arange(len(pdist)),np.array(pdist)))
        elif isinstance(pdist, scipy.stats.rv_discrete):
            self.pdist = pdist
        else:
            raise ValueError("""pdist must be array-like, scipy.stats.rv_discrete
                             object or one of the valid probability
                             distributions as a str""")   
                             
        return
    
    def sample_counts(self, *args, **kwargs):
        return self.pdist.rvs(*args, **kwargs)
    
    def random_walk(*args, **kwargs):
        """
        Simulate tunneling given some imaging time and a fixed tunneling rate. 

        Parameters
        ----------
        *args : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        return