#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   y0_plane_streamwise_vars.py
@Time    :   2025/02/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Read in y=0 plane statistics and plot streamwise (spanwise-averaged) 
             variables.
'''

import os
import sys
import numpy              as     np
import pandas             as     pd

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.params      import Params
from   vista.directories import Directories
from   vista.statistic   import StatisticData

# =============================================================================

casefolder = '/home/wencan/temp/250120'

# =============================================================================

dirs   = Directories( casefolder )

params = Params( dirs.case_para_file )

grid   = GridData( dirs.grid )
grid.read_grid()
block_list, _ = grid.select_sliced_blockgrids( 'Y', 0.001 )

stat2d = StatisticData( dirs.stat_yslice )
stat2d.read_statistic( block_list, vars_in=['pp','p'] )
stat2d.match_grid( block_list, grid, stat_type='Y' )
stat2d.drop_ghost( block_list )

for bl in stat2d.bl_clean:
    bl.df['p_fluc'] = np.sqrt(bl.df['pp']-bl.df['p']*bl.df['p'])/params.p_ref

df_all = pd.concat( [bl.df for bl in stat2d.bl_clean] )
df_all = df_all.groupby('x').mean()

print(df_all)


