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
import numpy              as     np
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.timer        import timer
from   vista.probe        import ProbeData
from   vista.probe        import ProbeFile
from   vista.params       import Params
from   vista.directories  import Directories
from   vista.directories  import create_folder
from   vista.tools        import create_linear_interpolator
from   vista.plot_setting import set_plt_rcparams
from   vista.plot_setting import set_plt_style

set_plt_rcparams()

# =============================================================================

case_dir  = '/home/wencan/temp/241030/' 
plot_all  = False
plot_stat = True

# =============================================================================

dirs     = Directories( case_dir )
probes   = ProbeFile( dirs.set_prb )

# -- enter working directory and read relevant parameters

os.chdir( dirs.prb_dir )

prb_files = dirs.probes
n_data = len( prb_files )
print(f"We are in directory:{dirs.prb_dir}\n")
print(f"We have got {n_data:5d} probes data.\n")

params    = Params( dirs.case_para_file )
p_ref     = params.p_ref
h         = params.H
x_imp     = params.x_imp
delta_0   = params.delta_0
# the following parameters are normalized by delta_0
x_sep     = params.x_sep
x_att     = params.x_att
x_pfmax   = params.x_pfmax
prb_withT = params.prb_withT

# -- create output directories

create_folder( dirs.pp_pre_ridge )
os.chdir( dirs.pp_pre_ridge )

# -- start read pressure at the ridge

x_locs = list()
pres   = list()

if not os.path.exists( 'pressure_ridge.pkl' ):
    
    with timer("Reading pressure data"):
        
# ----- first get the index of probes signals at snapshot time points.
#       (no matter which probe, as long as it contains the same time points like others)

        prb_data = ProbeData( prb_files[0], withT=prb_withT )
        prb_data.cleandata( t_start=20.0 )
        index = prb_data.time_index( np.linspace(20.0, 61.0, 4101) )
        
# ----- read pressure data at the ridge probe by probe

        for i in params.prb_ridge_index:
            
            probe = probes.probes[i-1]
            xyz = probe.xyz
            
            prb_data = ProbeData( prb_files[i-1], withT=prb_withT )
            prb_data.cleandata( t_start=20.0 )
            prb_data.df = prb_data.df.iloc[index]
             
            x_locs.append( xyz[0] )
            pres.append( np.array(prb_data.df['p']) )
            
            print(f"Reading pressure data at probe {i} at x={xyz[0]:.2f}.")

        times = np.array( prb_data.df['time'] )

        # -- save data into a pickle file
        with open( 'pressure_ridge.pkl', 'wb' ) as f:
            
            pickle.dump( [times, x_locs, pres], f )


# -- plot instantaneous pressure distribution at the ridge

with open( 'pressure_ridge.pkl', 'rb' ) as f:
    times, x_locs, pres = pickle.load( f )

pres = np.array( pres )
n_locs, n_time = pres.shape
x_locs = (np.array(x_locs) - x_imp) / delta_0

print(f"there are {n_locs} probes at the ridge, and {n_time} time frames.")

if plot_all:
    
    os.chdir( create_folder('./all_pressures') )
    
    for i in range( len(times) ):
        
        fig, ax = plt.subplots()
        
        set_plt_style( case='wall_pressure', ax=ax, fig=fig )
        
        pre = pres[:,i] / p_ref
        
        ax.plot( x_locs, pre )
        ax.set_title( f"Pressure distribution at t={times[i]:.2f}s" )
        
        plt.savefig( f"pressure_{i:06d}.png" )
        plt.close()
    
if plot_stat:
    
    df = pickle.load( open(dirs.pp_wall_proj+'/streamwise_vars.pkl', 'rb') )
    
    fig, ax = plt.subplots( )
    
    set_plt_style( case='wall_pressure', ax=ax, fig=fig )
    
    for i in range( 0, len(times), 3 ):
        
        ax.plot( x_locs[:-2], pres[:-2,i]/p_ref, alpha=0.01, color='gray' )
    
    ax.plot( df['x'], df['Cp'], label='avg', color='red' )
    
    f = create_linear_interpolator( df['x'], df['Cp'] )
    
    ax.plot([x_sep, x_att],[f(x_sep), f(x_att)], color='blue', marker='p',markersize=10,ls='')
    ax.plot([x_pfmax],[f(x_pfmax)], color='red', marker='p',markersize=10,ls='')

    plt.savefig( f"pressure_stat_every5.png" )
    plt.show()   
    plt.close()
    
        
        