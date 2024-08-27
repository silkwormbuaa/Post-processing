#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   snap_isosurface.py
@Time    :   2024/08/26
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   visualise the isosurface of shock and vortices
'''


import os
import sys
import time

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

import pyvista           as     pv

from   scipy.interpolate import griddata

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.plane_analy import save_sonic_line
from   vista.plane_analy import save_separation_line
from   vista.plane_analy import shift_coordinates

from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_slicez_stat

from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================
# option 
# =============================================================================

loc       = 0.0
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
    
#print(g.px, g.py, g.pz)

bbox = [-50.00, 125.0, 0.0, 50.0, -11.0, 11.0]
block_list = grid3d.select_blockgrids( bbox, mode='within' )

# - read in 3D snapshot file

with timer("read in 3d snapshot "):

    snap3d = Snapshot( snapshotfile )

    snap3d.verbose = False
    
    snap3d.grid3d = grid3d

    snap3d.read_snapshot( block_list )
    
#    snap3d.compute_gradients( block_list, grads, grid3d, buff=3 )

pvgrids = []

bl_df = snap3d.snap_data[0][5]

g = grid3d.g[snap3d.snap_data[0][0]-1]

grid = pv.RectilinearGrid(g.px,g.py,g.pz)

grid.cell_data['value'] = bl_df['u']

grid.plot()
#contours = grid.contour( 509 )
#
#pl = pv.Plotter()
#pl.add_mesh(grid)
#pl.add_mesh(contours)
#
#pl.show()

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()       