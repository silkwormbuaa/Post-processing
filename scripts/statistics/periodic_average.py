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
from   vista.directories import Directories
from   vista.params      import Params
from   vista.tools       import define_wall_shape
from   vista.plane_analy import pv_interpolate

case_folder = '/home/wencan/temp/241030/'
vars_read  = ['u','v','w','T']
bbox       = [-100,0,-2,10,-20,20]

dirs        = Directories( case_folder )

params      = Params( dirs.case_para_file )

grid       = GridData( dirs.grid )
grid.read_grid()

blocklist, _ = grid.select_sliced_blockgrids('X', -27.496,bbox=bbox)

stat_file    = dirs.sup_dir + '/stat_xslice_upstream.bin'
stat         = StatisticData(stat_file)
stat.grid3d  = grid
stat.read_statistic(blocklist,vars_in=vars_read)
stat.match_grid(blocklist, grid, add_to_df=True )
stat.compute_vars(blocklist,['mach'])

# periodic averaging

stat.spanwise_periodic_average( blocklist, ['u','v','w','mach'], params.period )

# interpolate into a whole cartesian grid

pz = np.linspace( 0.0, 2.6, 101, endpoint=True)
py = np.linspace(-0.6, 2.4, 201, endpoint=True)
px = np.array([0.0])

dataset = pv.MultiBlock( stat.create_vtk_multiblock(blocklist,['u','v','w','mach']) )
df      = pv_interpolate( dataset, ['u','v','w','mach'], [px,py,pz] )

v    = np.array( df['v']    ).reshape( (len(py),len(pz)) )
w    = np.array( df['w']    ).reshape( (len(py),len(pz)) )
u    = np.array( df['u']    ).reshape( (len(py),len(pz)) )
mach = np.array( df['mach'] ).reshape( (len(py),len(pz)) )
z    = np.array( df['z']    ).reshape( (len(py),len(pz)) )
y    = np.array( df['y']    ).reshape( (len(py),len(pz)) )

# visualization

fig, ax = plt.subplots(1,1,figsize=(8,6))

cs     = ax.contourf(z, y, v, levels=51, cmap='RdBu_r')
csnoic = ax.contour( z, y, mach, levels=[1.0], colors='lime', linewidths=2.0)

ax.streamplot(z,y, w, v, color='black', linewidth=0.5, density=4)

ywall = define_wall_shape(pz, casecode='241030', write=False)

ax.fill_between( pz, -0.6, y2=ywall, color='gray', zorder=10 )

ax.set_aspect('equal')

plt.show()