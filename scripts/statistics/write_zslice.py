#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   write_zslice.py
@Time    :   2024/10/28 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Read in a 3D statistics.bin file and write out mid-plane z-slice with
             all the statistics variables.
'''


import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.statistic   import StatisticData
from   vista.directories import Directories

# =============================================================================

casedir   = '/home/wencan/temp/smooth_mid/'
slic_type = 'Z'
loc       = 0.001

# =============================================================================

os.chdir( casedir)
dirs = Directories( casedir )

grid = GridData( dirs.grid )
grid.read_grid()
blocklist, index = grid.select_blockgrids( slic_type, loc )

stat3d = StatisticData( dirs.statistics )
stat3d.read_statistic( blocklist, vars_in=stat3d.full_vars )
stat3d.grid3d = grid

stat2d = stat3d.get_slice( slic_type, loc )
stat2d.write_statistic( dirs.stat_zslice )

