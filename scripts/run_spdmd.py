# -*- coding: utf-8 -*-
'''
@File    :   vista_spdmd.py
@Time    :   2023/05/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   sparsity promoting dmd script
             two parameters of gamma, rho should be provided.
             'python3 run_spdmd.py 10000. 10000.'
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


# =============================================================================
# Take gamma, rho from command line
# =============================================================================

arguments = sys.argv

if len( arguments ) == 1:
    
    print("Default gamma=10000. and rho=10000 will be used.\n")
    gamma = None
    rho = None

elif len( arguments ) == 3:
    
    gamma = float( arguments[1] )
    rho   = float( arguments[2] )
    print(f"Got {len(arguments)-1} parameters: gamma={gamma}",end='')
    print(f" rho={rho}\n")

else: raise ValueError("gamma and rho should be provided!")


# =============================================================================
# Computation
# =============================================================================

# Directory that stores Psq file

Lsep = 58.9713176

velocity = 507.0

snap_dir = os.getcwd()


# Set the instance and read in Pqs

paradmd = ParaDmd( snap_dir )

with timer('Read in Pqs file '):
    
    paradmd.read_Pqs()
    

# Compute spdmd

with timer('Compute sparsity promoting dmd '):
    
    
    paradmd.compute_spdmd( gamma, rho )
    
    print(f"The performance loss of SPDMD is : {paradmd.Ploss:8.3f} %.")
    

St = paradmd.freq * Lsep / velocity

paradmd.save_spdmd_result()

paradmd.save_ind_spmode()


# =============================================================================
# plot eigen values and amplitudes
# =============================================================================

''' 
numpy support using tuple as index to look for the elements.
e.g.: 
a. (2,1) will look for a 2-D array's element array[2,1].
b. ([2,3],[1]) will look for the elements array[2,1] and array[3,1].
c. ([1,2,3]) will look for the elements in a 1-D array arr[1:3].
'''

plot_eigens( paradmd.mu, 
             paradmd.mu[ tuple([paradmd.ind_spmode]) ],
             filename='eigens-.png')

plot_amp_st( St, 
             paradmd.alphas,
             amp2 = np.abs(paradmd.alphas_pol),
             filename='amplitude.png',
             hidesmall=False)

plot_amp_st( St, 
             paradmd.alphas,
             amp2 = np.abs(paradmd.alphas_pol),
             filename='amplitude-hidesmall-.png')
