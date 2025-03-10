#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   convert_2_vtk.py
@Time    :   2024/08/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   convert snapshot to vtk or szplt format
'''

import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.directories import create_folder

# =============================================================================

output_path = '/home/wencan/test/241030'
snap_file = '/home/wencan/temp/241030/snapshots/snapshot_01200380/snapshot.bin'
grid_file = '/home/wencan/temp/241030/results/inca_grid.bin'
output_filename = 'snapshot_01200380'
box  = [-30,999,-999.0,31.0,-999,999]     # box should also within the snapshot's range
vars_in = ['u','v','w','p','T']
vars_out =  ['u','v','w','p','T','Q_cr','vorticity','grad_rho_mod','div']
vtk   = True
szplt = False

# =============================================================================

# - read in grid data

G = GridData( grid_file )
G.read_grid()
blocks_list = G.select_blockgrids( box, mode='within' )

# - read in snapshot

with timer("read in snapshot"):
    
    snapshot = Snapshot( snap_file )
    snapshot.grid3d = G
    snapshot.verbose = False
    snapshot.read_snapshot( block_list=blocks_list, var_read=vars_in )
    
    snapshot.compute_gradients( block_list=blocks_list, grads=['grad_rho','vorticity','grad_rho_mod','div','Q_cr'] )

# - output

create_folder( output_path ); os.chdir( output_path )

if vtk:
    output_filename += '.vtm'
    with timer("write vtm"):
        snapshot.write_vtm( output_filename, vars=vars_out, block_list=blocks_list )

if szplt:
    output_filename += '.szplt'
    with timer("write szplt"):
        snapshot.write_szplt( output_filename, vars=vars_out, block_list=blocks_list )