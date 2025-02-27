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
import numpy             as     np
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.statistic   import StatisticData
from   vista.plane_analy import pv_interpolate


grid_file  = '/home/wencanwu/test/231124/results/inca_grid.bin'
stat_file  = '/home/wencanwu/test/231124/supplements/stat_xslice_000.bin'
vars_read  = ['u','v','w']
bbox       = [-100,0,-2,10,-20,20]

grid       = GridData(grid_file)
grid.read_grid()

blocklist, _ = grid.select_sliced_blockgrids('X', -27.496,bbox=bbox)

stat         = StatisticData(stat_file)
stat.grid3d  = grid
stat.read_statistic(blocklist,vars_in=vars_read)
stat.match_grid(blocklist, grid, add_to_df=True )

# periodic averaging

stat.spanwise_periodic_average( blocklist, ['u','v','w'], 1.3 )

# interpolate into a whole cartesian grid

pz = np.linspace(-2.6, 2.6, 101, endpoint=True)
py = np.linspace(-0.6, 2.4, 201, endpoint=True)
px = np.array([0.0])

dataset = pv.MultiBlock( stat.create_vtk_multiblock(blocklist,vars_read) )

p = pv.Plotter()
p.add_mesh( dataset, scalars='v')
p.show()

df      = pv_interpolate( dataset, vars_read, [px,py,pz] )

v = np.array( df['v'] ).reshape( (len(py),len(pz)) )
w = np.array( df['w'] ).reshape( (len(py),len(pz)) )
u = np.array( df['u'] ).reshape( (len(py),len(pz)) )
z = np.array( df['z'] ).reshape( (len(py),len(pz)) )
y = np.array( df['y'] ).reshape( (len(py),len(pz)) )

# visualization

fig, ax = plt.subplots(1,1,figsize=(8,6))

cs = ax.contourf(z,y,u, levels=100, cmap='jet')

ax.streamplot(z,y, w, v, color='black', linewidth=0.5, density=4)

plt.show()