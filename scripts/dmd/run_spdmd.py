#!/usr/bin/env python3
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
import time

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

import numpy             as     np

from   vista.timer       import timer

from   vista.paradmd     import ParaDmd

from   vista.colors      import colors  as col

from   vista.plot_style  import plot_eigens
from   vista.plot_style  import plot_amp_st
from   vista.plot_style  import plot_psi_st

from   vista.tools       import read_case_parameter

from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================
# Take gamma, rho from command line | read case parameters
# =============================================================================

arguments = sys.argv

if len( arguments ) == 1:
    
    print( "\nDefault", col.bg.green, col.fg.red, end='')
    print("gamma=10000, rho=10000 ", col.reset, "will be used.\n")
    
    gamma = None
    rho = None

elif len( arguments ) == 3 or len( arguments ) > 3 : # in case >run.out
    
    gamma = float( arguments[1] )
    rho   = float( arguments[2] )
    print(f"\nGot {len(arguments)-1} parameters: ", col.bg.green, end='')
    print(col.fg.red, f"gamma={gamma}, rho={rho}",  col.reset,"\n")

else: raise ValueError("gamma and rho should be provided!")

case_parameters = read_case_parameter( 'case_parameters' )

sys.stdout.flush()

# =============================================================================
# Computation
# =============================================================================

# Directory that stores Psq file

Lsep = float( case_parameters.get('Lsep') )

velocity = float( case_parameters.get('u_ref') )

snap_dir = os.getcwd()


# check if paradmd.pkl exists

if not os.path.exists( 'spdmd_result.pkl' ):

    # Set the instance and read in Pqs

    paradmd = ParaDmd( snap_dir )

    with timer('Read in Pqs file '):
        
        paradmd.read_Pqs()
        

    # Compute spdmd

    with timer('Compute sparsity promoting dmd '):
        
        
        paradmd.compute_spdmd( gamma, rho )
        
        print(f"\nThe performance loss of SPDMD is :",col.fg.green, end='')
        print(f" {paradmd.Ploss:8.3f} %.", col.reset, "\n")
        sys.stdout.flush()

    paradmd.St = paradmd.freq * Lsep / velocity

    paradmd.save_spdmd_result()

    paradmd.save_ind_spmode()

    sys.stdout.flush()
    
else: 
    
    paradmd = ParaDmd( snap_dir )
    
    paradmd.read_spdmd_result()


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

fmt = '.pdf'

eigen_dir = snap_dir + '/eigen_plots'

if not os.path.exists( eigen_dir ):
    os.mkdir( eigen_dir )
    print(f"Created directory {eigen_dir}.\n")
    
os.chdir( eigen_dir )

plot_eigens( paradmd.mu, 
             paradmd.mu[ tuple([paradmd.ind_spmode]) ],
             filename='a-eigens-' + fmt)

plot_amp_st( paradmd.St, 
             np.abs(paradmd.alphas),
             amp2 = np.abs(paradmd.alphas_pol),
             filename='a-amplitude'+fmt,
             hidesmall=False)

plot_amp_st( paradmd.St, 
             np.abs(paradmd.alphas),
             amp2 = np.abs(paradmd.alphas_pol),
             filename='a-amplitude-hidesmall'+fmt)

plot_psi_st( paradmd.St,
             np.abs(paradmd.psi),
             psi2 = np.abs(paradmd.psi_pol),
             filename='a-psi'+fmt)

plot_psi_st( paradmd.St,
             np.abs(paradmd.psi),
             psi2 = np.abs(paradmd.psi_pol),
             filename='a-psi-log'+fmt,
             xlim=[0.001,6.0])

growth = paradmd.beta
gmax = np.max(growth)
gmin = np.min(growth)

gray_scale = (growth - gmin)/( gmax - gmin )

plot_psi_st( paradmd.St,
             np.abs(paradmd.psi),
             psi2 = np.abs(paradmd.psi_pol),
             filename='a-psi-log-gray'+fmt,
             xlim=[0.001,6.0],
             gray=gray_scale)

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()      