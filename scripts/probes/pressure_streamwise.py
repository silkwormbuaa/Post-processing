#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   pressure_streamwise.py
@Time    :   2024/09/06 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Analyse the instantaneous streamwise pressure distribution evolution 
'''

import os
import sys
import pickle
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.timer       import timer
from   vista.probe       import ProbeData
from   vista.probe       import ProbeFile
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.tools       import read_case_parameter

# =============================================================================

case_dir = '/media/wencanwu/Seagate Expansion Drive1/temp/220927'

# =============================================================================

dirs     = Directories( case_dir )
probes   = ProbeFile( dirs.set_prb )

# -- enter working directory and read relevant parameters

os.chdir( dirs.prb_dir )

prb_files = dirs.probes
n_data = len( prb_files )
print(f"We are in directory:{dirs.prb_dir}\n")
print(f"We have got {n_data:5d} probes data.\n")

parameters = read_case_parameter( dirs.case_para_file )
p_ref   = float( parameters.get('p_ref') )
h       = float( parameters.get('H') )
x_imp   = float( parameters.get('x_imp') )
delta_0 = float( parameters.get('delta_0') )
prb_withT = True if parameters.get('prb_withT').lower() == 'true' else False

# -- create output directories

create_folder( dirs.pp_pre_ridge )
os.chdir( dirs.pp_pre_ridge )

# -- start read pressure at the ridge

x_locs = list()
pres   = list()

if not os.path.exists( 'pressure_ridge.pkl' ):
    
    with timer("Reading pressure data"):
        
        for i in range( len(probes.probes) ):
            
            probe = probes.probes[i]
            xyz = probe.xyz
            
            # find the probes at the ridges
            if not(abs(xyz[1]) < 0.001 and abs(xyz[2]) < 0.001):
                continue
            
            prb_data = ProbeData( prb_files[i], withT=prb_withT, step=40 )
            prb_data.cleandata( t_start=20.0 )
            
            x_locs.append( xyz[0] )
            pres.append( np.array(prb_data.df['p']) )
            
            print(f"Reading pressure data at probe {i} at x={xyz[0]:.2f}.")

        times = np.array( prb_data.df['time'] )

        # -- save data into a pickle file
        with open( 'pressure_ridge.pkl', 'wb' ) as f:
            
            pickle.dump( [times, x_locs, pres], f )
        
    
with timer("plot data"):
    
    with open( 'pressure_ridge.pkl', 'rb' ) as f:
        times, x_locs, pres = pickle.load( f )
    
    pres = np.array( pres )
    n_locs, n_time = pres.shape
    x_locs = (np.array(x_locs) - x_imp) / delta_0
    
    print(f"there are {n_locs} probes at the ridge, and {n_time} time frames.")
    
    print(x_locs[0], x_locs[-1])

    for i in range( len(times) ):
        
        fig, ax = plt.subplots( figsize=(15, 8) )
        
        pre = pres[:,i] / p_ref
        
        ax.plot( x_locs, pre )
        ax.set_title( f"Pressure distribution at t={times[i]:.2f}s" )
        ax.set_xlim(-12.5,12)
        ax.set_ylim(0.8, 2.5)
        
        plt.savefig( f"pressure_{i:06d}.png" )
        plt.close()
    
    
    
        
        