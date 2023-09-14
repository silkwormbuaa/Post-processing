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

datapath = "/home/wencanwu/my_simulation/temp/220927_lowRe/snapshots/snapshot_00699936"

# - read in grid file

gridfile = datapath.split('/snapshots')[0] + '/results/inca_grid.bin'

grid3d = GridData( gridfile )
grid3d.read_grid()

block_list, indx_slic = grid3d.select_sliced_blockgrids( slic_type, loc )

# - read in 3D snapshot file

with timer("read in 3d snapshot "):

    snapshotfile = datapath + '/snapshot.bin'

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
    
    x_slice = np.array( df_slice['x'] )
    y_slice = np.array( df_slice['y'] )
    u_slice = np.array( df_slice['u'] )
    
    x = np.linspace( -30.0, 100, 391 )
    y = np.linspace( 0.0, 40.0, 161 )
    
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
   
    
    





