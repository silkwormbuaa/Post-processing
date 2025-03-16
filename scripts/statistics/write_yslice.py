#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   write_zslice.py
@Time    :   2025/02/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Read in a 3D statistics.bin file and write out y-normal plane with
             all the statistics variables.
'''


import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder

# =============================================================================
# check memory estimation compared to the full statistics.bin file

casedir        = '/home/wencan/temp/241030/'
slic_type      = 'Y'
default_output = False

loc_specified  = 0.1
out_specified  = 'y_0.1/stat_y_0.1.bin'


# =============================================================================

os.chdir( casedir )
dirs = Directories( casedir )

if default_output:
    loc         = 0.001
    output_file = dirs.stat_yslice
else:
    loc         = loc_specified
    output_file = dirs.pp_statistics + '/' + out_specified
    os.chdir( dirs.pp_statistics )
    create_folder( out_specified.split('/')[0] )
    
grid = GridData( dirs.grid )
grid.read_grid()
blocklist, index = grid.select_sliced_blockgrids( slic_type, loc )

stat3d = StatisticData( dirs.statistics )
stat3d.read_statistic( blocklist, vars_in=stat3d.full_vars )
stat3d.grid3d = grid

stat2d = stat3d.get_slice( slic_type, loc )
stat2d.write_statistic( output_file )

