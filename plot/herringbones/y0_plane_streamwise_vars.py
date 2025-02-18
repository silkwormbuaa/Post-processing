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
import pickle
import numpy              as     np
import pandas             as     pd
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.line        import LineData
from   vista.params      import Params
from   vista.directories import Directories
from   vista.statistic   import StatisticData
# =============================================================================

casefolder = '/home/wencan/temp/250120'

file0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/wall_projection/streamwise_vars.pkl'

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

df_all = df_all.groupby('x').mean().reset_index() # reset_index() is necessary to make x a column

print(df_all)

# =============================================================================

line = LineData()
with open( file0, 'rb' ) as f:    line.df = pickle.load( f )

# =============================================================================

plt.rcParams['font.size'] = 20
fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.plot((df_all['x']-50.4)/5.2, df_all['p_fluc'], label='p_fluc', lw=2,color='black')
ax.plot( line.df['x'], line.df['p_fluc'], label='smoothwall', color='r', lw=2 )

ax.set_xlabel(r'$(x-x_\delta)/\delta_0$')
ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_\infty$")

ax.set_xlim([-20,10])

plt.show()


