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
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.statistic   import StatisticData
from   vista.directories import Directories

# =============================================================================
# check memory estimation compared to the full statistics.bin file

casedir        = '/home/wencan/temp/241030/'
slic_type      = 'Y'

locs           = [0.025]
locs           = np.array(locs) * 5.2
outfiles       = [f'/stat_yslice_{i:03d}.bin' for i in range(len(locs))]

# =============================================================================

os.chdir( casedir )
dirs = Directories( casedir )

grid = GridData( dirs.grid )
grid.read_grid()

for i, loc in enumerate( locs ):
    
    output_file = dirs.sup_dir + outfiles[i]

    blocklist, index = grid.select_sliced_blockgrids( slic_type, loc )

    stat3d = StatisticData( dirs.statistics )
    stat3d.read_statistic( blocklist, vars_in=stat3d.full_vars )
    stat3d.grid3d = grid

    stat2d = stat3d.get_slice( slic_type, loc )
    stat2d.write_statistic( output_file )

