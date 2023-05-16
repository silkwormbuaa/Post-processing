# -*- coding: utf-8 -*-
'''
@File    :   vista_spdmd.py
@Time    :   2023/05/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   sparsity promoting dmd script
'''

import os

import sys

import pickle

import numpy             as     np

import pandas            as     pd

sys.path.append('..')

from   utils.timer       import timer

from   vista.paradmd     import ParaDmd

from   utils.plot_style  import plot_eigens

from   utils.plot_style  import plot_amp_st

# Directory that stores Psq file
snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots'

#snap_dir = '/home/wencanwu/my_simulation/temp/220926_lowRe/snapshots/Z_slice_results'
# gammas that are going to be sweep over

gammas = np.arange(10000000,11000000,1000000)


# Set the instance and read in Pqs

paradmd = ParaDmd( snap_dir )

with timer('Read in Pqs file '):
    
    paradmd.read_Pqs()
    

# Compute spdmd

with timer('Compute sparsity promoting dmd '):
    
    paradmd.compute_spdmd( gammas )
    



Lsep = 50.069

velocity = 507.0

St = paradmd.freq * Lsep / velocity

for i in range(len(gammas)):

    plot_eigens( paradmd.mu, paradmd.mu[paradmd.ind_sel[i]] )
    
    plot_amp_st( St, paradmd.alphas ,
                amp2 = np.abs(paradmd.alphas_pol[i,:]) )