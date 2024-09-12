#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   snap_wall_vars.py
@Time    :   2024/09/08 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   using pyvista to extract variables at the wall
'''

import os
import gc
import sys
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot

# =============================================================================

snap_file = '/media/wencan/Expansion/temp/220927/supplements/wall_dist/snapshot_01329504/snapshot.bin'
grid_file = '/media/wencan/Expansion/temp/220927/results/inca_grid.bin'
bbox      = [-30,999,-3,31.0, 0,999]
vars_in   = ['u','p','wd']
# =============================================================================

grid3d = GridData( grid_file )
grid3d.read_grid()
block_list = grid3d.select_blockgrids( bbox, mode='within' )

snap3d = Snapshot( snap_file )
snap3d.grid3d = grid3d
# snap3d.verbose = True
snap3d.read_snapshot( block_list=block_list, var_read=vars_in )

dataset = pv.MultiBlock( snap3d.create_vtk_multiblock(vars=vars_in,block_list=block_list,buff=2))

point_data = dataset.cell_data_to_point_data().combine()
point_data.set_active_scalars('wd')

wallsurface = point_data.contour( [0.01] )
wallsurface.set_active_scalars('p')

p = pv.Plotter()
cmap = plt.get_cmap('coolwarm',51)

p.add_mesh( wallsurface, cmap=cmap, show_scalar_bar=True)

p.add_axes()

p.show()

p.close()


