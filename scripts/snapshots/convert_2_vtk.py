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

output_path = 'xxxx'
snap_file = 'xxxx/snapshot.bin'
grid_file = 'xxxx/grid.bin'
output_filename = 'snapshot_0000.vtm'
box  = [-999,999,-999,999,-999,999]
vars = ['u']

# - read in grid data

G = GridData( grid_file )
G.read_grid()

# - read in snapshot

with timer("read in snapshot"):
    
    snapshot = Snapshot( snap_file )
    snapshot.grid3d = G
    snapshot.verbose = False

    snapshot.read_snapshot( var_read=vars )

# - output

with timer("write vtm"):
    create_folder( output_path ); os.chdir( output_path )
    snapshot.write_vtm( output_filename, vars )