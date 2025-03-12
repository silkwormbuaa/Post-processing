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
import pyvista            as     pv
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.line        import LineData
from   vista.params      import Params
from   vista.directories import Directories
from   vista.statistic   import StatisticData
# =============================================================================

cases = ['250120','250218','250304']
colors = ['r','b','g']

file0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/wall_projection/streamwise_vars.pkl'

# =============================================================================

plt.rcParams['font.size'] = 20
fig, ax = plt.subplots(1,1,figsize=(8,6))

for i, case in enumerate(cases):

    casefolder = '/home/wencan/temp/' + case

    dirs   = Directories( casefolder )

    params = Params( dirs.case_para_file )

    grid   = GridData( dirs.grid )
    grid.read_grid()
    block_list, _ = grid.select_sliced_blockgrids( 'Y', 0.001 )

    stat2d = StatisticData( dirs.stat_yslice )
    stat2d.grid3d = grid
    stat2d.read_statistic( block_list, vars_in=['pp','p'] )
    stat2d.match_grid( block_list, grid )
    stat2d.drop_ghost( block_list )

    dataset = pv.MultiBlock( stat2d.create_vtk_multiblock(block_list, vars=['pp','p']) )

    dataset = dataset.cell_data_to_point_data().combine()

    # - local pressure fluctuation

    for bl in stat2d.bl_clean:
        bl.df['p_fluc'] = np.sqrt(bl.df['pp']-bl.df['p']*bl.df['p'])/params.p_ref
        
    # - spanwise dispersive pressure fluctuation

    df_xz = pd.concat( [bl.df for bl in stat2d.bl_clean] )

    df_xz['p_fluc']    = np.sqrt(df_xz['pp']-df_xz['p']**2)/params.p_ref
    df_xz['p_zmean']   = df_xz.groupby('x')['p'].transform('mean')
    df_xz['p_fluc_total'] = np.sqrt(np.abs(df_xz['pp']-df_xz['p_zmean']**2))/params.p_ref
    df_xz['p_2']       = df_xz['p']**2
    df_xz['p_2_zmean'] = df_xz.groupby('x')['p_2'].transform('mean')

    df_xz['p_disp']    = np.sqrt( df_xz['p_2_zmean'] - df_xz['p_zmean']**2 )/params.p_ref

    df_zmean = df_xz.groupby('x').mean().reset_index() # reset_index() is necessary to make x a column

    # dataset['p_fluc']          = df_xz['p_fluc'].values
    # dataset['p_fluc_total']    = df_xz['p_fluc_total'].values
    # dataset['p_disp']          = df_xz['p_disp'].values
    # dataset['p_zmean']         = df_xz['p_zmean'].values
    # dataset['p_2_zmean']       = df_xz['p_2_zmean'].values
    # dataset['p_2']             = df_xz['p_2'].values

    # print( df_xz )

    # print(df_zmean)

    # p = pv.Plotter()

    # clipbox = [-100,0.0,-20,20,-20,20]

    # dataset = dataset.clip_box(clipbox, invert=False)

    # p.add_mesh( dataset, scalars='p_disp', cmap='RdBu_r' )
    # p.add_axes()
    # # p.show()
    # p.close()

# =============================================================================

    ax.plot((df_zmean['x']-50.4)/5.2, df_zmean['p_disp'], label='p_disp', lw=1, color=colors[i])
    
    #ax.plot((df_zmean['x']-50.4)/5.2, df_zmean['p_fluc'], label='p_disp', lw=1, color=colors[i])
    #ax.plot( line.df['x'], line.df['p_fluc'], label='smoothwall', color='r', lw=2 )

# =============================================================================

# line = LineData()
# with open( file0, 'rb' ) as f:    line.df = pickle.load( f )

ax.set_xlabel(r'$(x-x_\delta)/\delta_0$')
ax.set_ylabel(r"$p_d$")

ax.set_xlim([-20,10])

plt.show()
plt.close()
