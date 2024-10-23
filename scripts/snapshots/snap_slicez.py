#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   snap_slicez.py
@Time    :   2023/09/14 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import time
import numpy             as     np
from   scipy.interpolate import griddata

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.snapshot    import Snapshot
from   vista.params      import Params
from   vista.plane_analy import save_sonic_line
from   vista.plane_analy import save_separation_line
from   vista.plane_analy import shift_coordinates
from   vista.plane_analy import compute_DS
from   vista.plot_style  import plot_slicez_stat
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )


# =============================================================================
# option 
# =============================================================================

slic_type = 'Z'
loc       = 0.0
grads     = ['grad_rho']

# =============================================================================
# read in grid file and snapshot file, then get slice dataframe
# =============================================================================

datapath = os.getcwd()
snapshotfile = datapath + '/snapshot.bin'
parametersfile = datapath.split('/snapshots')[0] + '/case_parameters'

# - read in grid file

gridfile = datapath.split('/snapshots')[0] + '/results/inca_grid.bin'

grid3d = GridData( gridfile )
grid3d.read_grid()

bbox = [-50.00, 125.0, -1.30, 50.0, -11.0, 11.0]
block_list, indx_slic = grid3d.select_sliced_blockgrids( slic_type, loc, bbox)

# - read in 3D snapshot file

with timer("read in 3d snapshot "):

    snap3d = Snapshot( snapshotfile )

    snap3d.verbose = False
    
    snap3d.grid3d = grid3d

    snap3d.read_snapshot( block_list )
    
    snap3d.compute_gradients( block_list, grads, buff=3 )
    

with timer("get slice dataframe "):
    
    df_slice = snap3d.get_slice_df( slic_type, loc )

# =============================================================================
# Interpolate and plot
# =============================================================================

with timer("Interpolate and plot "):
    
    params  = Params( parametersfile )
    delta   = params.delta_0
    h_ridge = params.H
    h_md    = params.H_md
    x_imp   = params.x_imp
    
    df_slice = shift_coordinates( df_slice, delta, h_ridge, h_md, x_imp )  
    
    x_slice = np.array( df_slice['xs'] )
    y_slice = np.array( df_slice['y_scale'] )
    
    u_slice = np.array( df_slice['u'] )
    grad_rho_slice = np.array( df_slice['grad_rho'] )
    
    x = np.linspace( -14, 12.5, 1801 )
    if loc == 0.0:
        y = np.linspace( 0.01, 8, 961)
    else:
        y = np.linspace( -0.1, 8, 325 )
    
    xx,yy = np.meshgrid( x, y )
    
    u = griddata( (x_slice,y_slice), u_slice,
                  (xx,yy), method='linear')

    grad_rho = griddata( (x_slice,y_slice), grad_rho_slice,
                         (xx,yy), method='linear')
    
    DS = compute_DS( grad_rho, min=0.0, max=2.0 )

    save_separation_line( xx,yy,u )
    
    cbar = r'$\nabla{\rho}$'
    cbar_levels = np.linspace( 0.0, 0.75*np.max(grad_rho),51)
    
    plot_slicez_stat( xx,yy,grad_rho,
                      filename='schlieren',
                      col_map='Greys',
                      cbar_label=cbar,
                      separation='separationlines.pkl',
                      sonic=False,
                      cbar_levels=cbar_levels)

    cbar = r'$DS$'
    cbar_levels = np.linspace( 0.0, 0.8,33)
    
    plot_slicez_stat( xx,yy,DS,
                      filename='DS_schlieren',
                      col_map='Greys_r',
                      cbar_label=cbar,
                      separation='separationlines.pkl',
                      sonic=False,
                      cbar_levels=cbar_levels,
                      x_lim=[-13,10],
                      y_lim=[0,8],
                      pure=False)
    
    cbar = r'$u/u_{\infty}$'
    cbar_levels = np.linspace( -0.2, 1, 37)
    cbar_ticks  = np.linspace( -0.2, 1, 7)
    
    plot_slicez_stat( xx,yy,u/507,
                      filename='streamwise_velocity',
                      col_map='coolwarm',
                      cbar_label=cbar,
                      separation='separationlines.pkl',
                      sonic=False,
                      cbar_levels=cbar_levels,
                      cbar_ticks=cbar_ticks,
                      x_lim=[-13,10],
                      y_lim=[0,8],
                      pure=False)

    
# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()       
    
   
    
    





