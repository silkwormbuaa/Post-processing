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
import pyvista           as     pv

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.timer       import timer
from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.plane_analy import save_sonic_line
from   vista.plane_analy import save_separation_line
from   vista.plane_analy import shift_coordinates
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================
# option 
# =============================================================================

loc       = 0.0
grads     = ['grad_rho']
bbox = [-50.00, 125.0, 0.0, 50.0, -11.0, 11.0]

# =============================================================================
# read in grid file and snapshot file, then get slice dataframe
# =============================================================================

datapath = os.getcwd()
snapshotfile = datapath + '/snapshot.bin'

# - read in grid file

gridfile = datapath.split('/snapshots')[0] + '/results/inca_grid.bin'
grid3d = GridData( gridfile )
grid3d.read_grid()
    
block_list = grid3d.select_blockgrids( bbox, mode='within' )

# - read in 3D snapshot file

with timer("read in 3d snapshot "):

    snap3d = Snapshot( snapshotfile )
    snap3d.verbose = False
    snap3d.grid3d = grid3d
    snap3d.read_snapshot( block_list )
    snap3d.compute_gradients( block_list, grads )

# - generate the vtk dataset

snap3d.create_vtk_multiblock( vars=grads, block_list=block_list )

# =============================================================================



print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()       