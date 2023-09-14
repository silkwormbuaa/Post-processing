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

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

import numpy             as     np

from   scipy.interpolate import griddata

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.plane_analy import save_sonic_line

from   vista.plane_analy import save_separation_line

from   vista.plane_analy import shift_coordinates

from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_slicez_stat

import matplotlib.pyplot as plt

# =============================================================================
# option 
# =============================================================================

slic_type = 'Z'
loc       = 0.0

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

block_list, indx_slic = grid3d.select_sliced_blockgrids( slic_type, loc )

# - read in 3D snapshot file

with timer("read in 3d snapshot "):

    snap3d = Snapshot( snapshotfile )

    snap3d.verbose = False
    
    snap3d.grid3d = grid3d

    snap3d.read_snapshot( block_list )

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
    
    x = np.linspace( -18, 12, 301 )
    y = np.linspace( 0.0, 7.5, 161 )
    
    xx,yy = np.meshgrid( x, y )
    
    u = griddata( (x_slice,y_slice), u_slice,
                  (xx,yy), method='linear')

   
    save_separation_line( xx,yy,u )
    
    cbar = 'u'
    plot_slicez_stat( xx,yy,u,
                      filename='utest',
                      cbar_label=cbar,
                      separation=True,
                      sonic=False)
   
    
    





