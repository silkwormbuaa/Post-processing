#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   extract_surface.py
@Time    :   2025/02/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Extract a surface from statistics data and then perform some integrations.
'''

import os
import sys
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.tools       import get_filelist

# =============================================================================

casedir = '/home/wencan/temp/250120/'
vars_in = ['rho','mu','p']
bbox    = [ -57.92, -48.70, -3.0, 1.0, -11, 11]

# =============================================================================

os.chdir( casedir )
dirs = Directories( casedir )

grid = GridData( dirs.grid )
grid.read_grid()

block_list = grid.select_blockgrids( bbox, mode='within' )

wd_file = get_filelist( dirs.wall_dist, 'snapshot.bin' )[0]
wd_snap = Snapshot( wd_file )
wd_snap.read_snapshot( block_list, var_read=['wd'] )

stat3d = StatisticData( dirs.statistics )
stat3d.read_statistic( block_list, vars_in=vars_in )
stat3d.grid3d = grid
stat3d.match_grid(block_list, grid)

for num in block_list:
    stat3d.bl[num-1].df['wd'] = wd_snap.snap_data[num-1].df['wd']

# =============================================================================

dataset = pv.MultiBlock(stat3d.create_vtk_multiblock( block_list, vars_in + ['wd'] ))

point_data = dataset.cell_data_to_point_data().combine()
point_data.set_active_scalars('wd')

# - extract surface

surface = point_data.contour( [00.0] )
surface.set_active_scalars( 'rho' )

integrated_data = surface.integrate_data()

print(integrated_data.array_names)

print(f"integrated rho over area is {integrated_data['rho'][0]}." )
print(f"integrated mu  over area is {integrated_data['mu'][0]}." )
print(f"integrated area is          {integrated_data['Area'][0]}.")
print(f"average rho is              {integrated_data['rho'][0]/integrated_data['Area'][0]}.")
print(f"average mu  is              {integrated_data['mu'][0]/integrated_data['Area'][0]}.")

p = pv.Plotter()
cmp = plt.get_cmap('coolwarm',51)

p.add_mesh( surface, scalars='p', cmap=cmp, show_scalar_bar=True )
p.add_axes()
p.show()
