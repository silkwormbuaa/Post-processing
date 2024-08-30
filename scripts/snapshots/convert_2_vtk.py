#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   convert_2_vtk.py
@Time    :   2024/08/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.directories import create_folder

output_path = '/home/wencanwu/test/231124'
snap_file = '/media/wencanwu/Seagate Expansion Drive1/temp/231124/snapshots/snapshot_02920790/snapshot.bin'
grid_file = '/media/wencanwu/Seagate Expansion Drive1/temp/231124/results/inca_grid.bin'
output_filename = 'snapshot_02920790.vtm'
box  = [-30,999,-999.0,31.0,0.0,999]     # box should also within the snapshot's range
vars = ['u','v','w','p','T']

# - read in grid data

G = GridData( grid_file )
G.read_grid()
blocks_list = G.select_blockgrids( box, mode='within' )

# - read in snapshot

with timer("read in snapshot"):
    
    snapshot = Snapshot( snap_file )
    snapshot.grid3d = G
    snapshot.verbose = False
    snapshot.read_snapshot( block_list=blocks_list, var_read=vars )
    
    snapshot.compute_gradients( block_list=blocks_list, grads=['grad_rho','vorticity','grad_rho_mod','div'] )

# - output

with timer("write vtm"):
    create_folder( output_path ); os.chdir( output_path )
    snapshot.write_vtm( output_filename, ['u','p','grad_rho_mod','div'], block_list=blocks_list )
