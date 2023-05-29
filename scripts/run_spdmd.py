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

source_dir = os.path.dirname( os.path.dirname( os.path.realpath(__file__) ))
sys.path.append( source_dir )

import numpy             as     np

from   vista.timer       import timer

from   vista.paradmd     import ParaDmd

from   vista.plot_style  import plot_eigens

from   vista.plot_style  import plot_amp_st

# Directory that stores Psq file
#snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots'

snap_dir = os.getcwd()
# gammas that are going to be sweep over

# gammas = np.arange(81500,801600,100)
gammas = np.array([ 450000000. ])

# Set the instance and read in Pqs

paradmd = ParaDmd( snap_dir )

with timer('Read in Pqs file '):
    
    paradmd.read_Pqs()
    

# Compute spdmd

with timer('Compute sparsity promoting dmd '):
    
    paradmd.compute_spdmd( gammas, rho=10000. )
    


Lsep = 50.069

velocity = 507.0

St = paradmd.freq * Lsep / velocity

for i in range(len(gammas)):

    plot_eigens( paradmd.mu, paradmd.mu[paradmd.ind_sel[i]] )
    
    plot_amp_st( St, paradmd.alphas ,
                amp2 = np.abs(paradmd.alphas_pol[i,:]) )
