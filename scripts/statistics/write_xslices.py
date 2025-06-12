#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   write_zslice.py
@Time    :   2024/11/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Read in a 3D statistics.bin file and write out a few x-slices.
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

casedir   = '/home/wencan/temp/smooth_mid/'
slic_type = 'X'
locs      = [-14.98,-8.53,-7.33,1.89,8.0]
locs      = np.array(locs) * 5.2 + 50.4
outfiles  = [f'/stat_xslice_{i:03d}.bin' for i in range(len(locs))]

#locs      = [-20.0]
#outfiles  = ['/stat_xslice_upstream.bin']

# =============================================================================

os.chdir( casedir)
dirs = Directories( casedir )

grid = GridData( dirs.grid )
grid.read_grid()

for i, loc in enumerate( locs ):
    
    outfile = dirs.sup_dir + outfiles[i]

    blocklist, index = grid.select_sliced_blockgrids( slic_type, loc )

    stat3d = StatisticData( dirs.statistics )
    stat3d.read_statistic( blocklist, vars_in=stat3d.full_vars )
    stat3d.grid3d = grid

    stat2d = stat3d.get_slice( slic_type, loc )
    stat2d.write_statistic( outfile )

