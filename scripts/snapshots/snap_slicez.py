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
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   scipy.interpolate import griddata

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.plane_analy import save_sonic_line
from   vista.plane_analy import save_separation_line
from   vista.plane_analy import shift_coordinates

from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_slicez_stat


# =============================================================================
# option 
# =============================================================================

slic_type = 'Z'
loc       = 0.65
grads     = ['schlieren']
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
    
    snap3d.compute_gradients( block_list, grads, grid3d, buff=3 )
    

with timer("get slice dataframe "):
    
    df_slice = snap3d.get_slice_df( slic_type, loc )

# =============================================================================
# Interpolate and plot
# =============================================================================

with timer("Interpolate and plot "):
    
    parameters = read_case_parameter( parametersfile )
    delta   = float( parameters.get('delta_0') )
    h_ridge = float( parameters.get('H') )
    h_md    = float( parameters.get('H_md') )
    x_imp   = float( parameters.get('x_imp') )
    
    df_slice = shift_coordinates( df_slice, delta, h_ridge, h_md, x_imp )  
    
    x_slice = np.array( df_slice['xs'] )
    y_slice = np.array( df_slice['y_scale'] )
    
    u_slice = np.array( df_slice['u'] )
    grad_rho_slice = np.array( df_slice['grad_rho'] )
    
    x = np.linspace( -17.5, 12.5, 601 )
    y = np.linspace( -0.1, 8, 325 )
    
    xx,yy = np.meshgrid( x, y )
    
    u = griddata( (x_slice,y_slice), u_slice,
                  (xx,yy), method='linear')

    grad_rho = griddata( (x_slice,y_slice), grad_rho_slice,
                         (xx,yy), method='linear')
   
    save_separation_line( xx,yy,u )
    
    cbar = r'$\nabla{\rho}$'
    cbar_levels = np.linspace( 0.0, 0.75*np.max(grad_rho),51)
    
    plot_slicez_stat( xx,yy,grad_rho,
                      filename='schlieren',
                      col_map='Greys',
                      cbar_label=cbar,
                      separation=True,
                      sonic=False,
                      cbar_levels=cbar_levels)
   
    
    





