#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   periodic_average.py
@Time    :   2025/02/25 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import pyvista           as     pv

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.statistic   import StatisticData
from   vista.grid        import GridData


grid_file  = '/home/wencanwu/test/231124/results/inca_grid.bin'
stat_file  = '/home/wencanwu/test/231124/supplements/stat_xslice_000.bin'
vars_read  = ['v']
bbox       = [-100,0,-2,10,-20,20]

grid       = GridData(grid_file)
grid.read_grid()

blocklist, _ = grid.select_sliced_blockgrids('X', -27.496,bbox=bbox)

stat         = StatisticData(stat_file)
stat.grid3d  = grid
stat.read_statistic(blocklist,vars_in=vars_read)
stat.match_grid(blocklist, grid, add_to_df=True )

# periodic averaging

stat.spanwise_periodic_average( blocklist, 'v', 1.3 )

#sys.exit()

# visualization

dataset = pv.MultiBlock( stat.create_vtk_multiblock(blocklist,vars_read) )

p = pv.Plotter()

dataset.set_active_scalars('v')

p.add_mesh(dataset)

p.show()

