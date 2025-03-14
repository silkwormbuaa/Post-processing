#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   spanwise_average.py
@Time    :   2025/03/14 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   just do direct ensemble average for spanwise direction,
             (no fluid/solid cells are considered).
             Results are save into a pickled pyvista dataset.
'''

import os
import sys
import pickle
import numpy             as     np
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.plane_analy import pv_interpolate

# =============================================================================

case_dir   = '/home/wencan/temp/241030/'
vars_read  = ['u','v','w','T']
bbox       = [-57.5,120,0,100,-20,20]

# =============================================================================

dirs = Directories( case_dir )

os.chdir( create_folder(dirs.pp_z_average) )

if not os.path.exists( dirs.pp_z_average + '/z_average.vtmb' ):

    print("z_average.pkl not found, start to compute spanwise average...\n")
    
    grid = GridData( dirs.grid ) 
    grid.read_grid()
    blocklist = grid.select_blockgrids( bbox=bbox, mode='within' )

    stat = StatisticData( dirs.statistics )
    stat.grid3d = grid
    stat.read_statistic( blocklist, vars_read )
    stat.compute_vars( blocklist, vars_new=['mach'] )
    dataset = pv.MultiBlock(stat.spanwise_average( blocklist, vars_read+['mach'] ))
    dataset.save( dirs.pp_z_average + '/z_average.vtmb', binary=True )

else:
    
    print("Found z_average.vtmb, loading it...\n")
    
    dataset = pv.read( dirs.pp_z_average + '/z_average.vtmb' )

# - post-process

px = np.linspace( -53.6, 100, 201, endpoint=True )
py = np.linspace(   0.0, 60,  121, endpoint=True )
pz = np.array([0.0])

dataset = dataset.cell_data_to_point_data().combine()

df = pv_interpolate( dataset, vars_read+['mach'],[px,py,pz])

x     = np.array( df['x']    ).reshape( (len(py),len(px)) )
y     = np.array( df['y']    ).reshape( (len(py),len(px)) )
v     = np.array( df['v']    ).reshape( (len(py),len(px)) )
w     = np.array( df['w']    ).reshape( (len(py),len(px)) )
u     = np.array( df['u']    ).reshape( (len(py),len(px)) )
mach  = np.array( df['mach'] ).reshape( (len(py),len(px)) )

fig, ax = plt.subplots(1,1,figsize=(8,6))
c    = ax.contourf( x, y, mach, levels=np.linspace(0.0,2.5,51), cmap='coolwarm', extend='both')
csep = ax.contour(  x, y, u,    levels=[0.0], colors='yellow', linewidths=1.5)
cson = ax.contour(  x, y, mach, levels=[1.0], colors='green',  linewidths=1.5)
ax.set_aspect('equal')

plt.show()
plt.close()

# def visualize( dataset, var ):
#     p = pv.Plotter()
#     p.add_mesh( dataset, scalars=var, show_edges=False, cmap='coolwarm')
#     p.add_axes()
#     p.view_vector([0.0,0.0,1.0], viewup=[0.0,1.0,0.0])
#     p.show()
#     p.close()

# visualize( dataset, 'mach' )

