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

from  vista.statistic    import StatisticData
from  vista.grid         import GridData
from  vista.timer        import timer
from  vista.directories  import create_folder

output_path = 'xxxx'
stat_file = 'xxxx/statistics.bin'
grid_file = 'xxxx/grid.bin'
box  = [-999,999,-999,999,-999,999]
vars = ['u']

# - read in grid data

G = GridData( grid_file )
G.read_grid()

# - determine which blocks to fill

block_list = G.select_blockgrids( box )

# - read in statistics

S = StatisticData( stat_file )

with timer("read in statistics"):
    with open( stat_file, 'rb' ) as f:
        S.verbose = True
        S.read_stat_header( f )
        S.read_stat_body( f, block_list, vars )

# - prepare grid 

S.match_grid( block_list, G, add_to_df=False )
S.grid3d = G

# - output

with timer("write vtm"):
    create_folder( output_path ); os.chdir( output_path )
    S.write_vtm( 'statistics.vtm', vars, block_list )