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
import pandas            as     pd
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

plt.rcParams["text.usetex"]         = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
plt.rcParams['font.family']         = "Times New Roman"
plt.rcParams['font.size']           = 30

# =============================================================================

case_folder = '/home/wencan/temp/220927/'
vars_read   = ['u','v','w','T']
bbox        = [-100,0,-2,10,-20,20]
streamline  = False
D_norm      = False

figname     = 'zoomin_v'
if streamline: figname += '_streamline'
if D_norm:     figname += '_D'
figname    += '.png'

# =============================================================================

dirs         = Directories( case_folder )
params       = Params( dirs.case_para_file )
grid         = GridData( dirs.grid )

grid.read_grid()
blocklist, _ = grid.select_sliced_blockgrids('X', -53.6, bbox=bbox)

stat_file    = dirs.sup_dir + '/stat_xslice_upstream.bin'
stat         = StatisticData(stat_file)
stat.grid3d  = grid
stat.read_statistic(blocklist,vars_in=vars_read)
stat.match_grid(blocklist, grid, add_to_df=True )
stat.compute_vars(blocklist,['mach'])

# periodic averaging

stat.spanwise_periodic_average( blocklist, ['u','v','w','mach'], params.D )

# interpolate into a whole cartesian grid

pz = np.linspace( 0.0, 1.3, 81, endpoint=True)
py = np.linspace(-0.6, 2.4, 201, endpoint=True)
px = np.array([0.0])

dataset = pv.MultiBlock( stat.create_vtk_multiblock(blocklist,['u','v','w','mach']) )
dataset = dataset.cell_data_to_point_data()
df      = pv_interpolate( dataset, ['u','v','w','mach'], [px,py,pz] )

v     = np.array( df['v']    ).reshape( (len(py),len(pz)) )/params.u_ref
w     = np.array( df['w']    ).reshape( (len(py),len(pz)) )/params.u_ref
u     = np.array( df['u']    ).reshape( (len(py),len(pz)) )/params.u_ref
mach  = np.array( df['mach'] ).reshape( (len(py),len(pz)) )

if D_norm: 
    length_unit = params.D
    y_bottom    = -0.48
    x_lim       = [0.0,1.0]
    y_lim       = [-0.48,1.6] 
    x_ticks     = [0.0,0.5,1.0]
    x_ticklabels= [r'$0.0$',r'$0.5$',r'$1.0$']
    y_ticks     = [-0.4,0.0,0.4,0.8,1.2,1.6]
    x_label     = r'$z/D$'
    y_label     = r'$y/D$'
    loc_tag     = [0.95, 1.4]
else     : 
    length_unit = params.delta_0
    y_bottom    = -0.12
    x_lim       = [0.0,0.25]
    y_lim       = [-0.12,0.4]
    x_ticks     = [0.0,0.125,0.25]
    x_ticklabels= [r'$0.00$',r'$0.125$',r'$0.25$']
    y_ticks     = [-0.1,0.0,0.1,0.2,0.3,0.4]
    x_label     = r'$z/\delta_0$'
    y_label     = r'$y/\delta_0$'
    loc_tag     = [0.23, 0.35]
 
ywall = define_wall_shape(pz, casecode=params.casecode, write=False)/length_unit
z     = np.array( df['z']    ).reshape( (len(py),len(pz)) )/length_unit
y     = np.array( df['y']    ).reshape( (len(py),len(pz)) )/length_unit

# visualization

fig = plt.figure(figsize=(8,6))
ax  = fig.add_axes([0.1,0.2,0.95,0.7])

cbar_levels = np.linspace( -2.0, 2.0, 51)
cbar_ticks  = np.linspace( -2.0, 2.0, 5)

cs     = ax.contourf(z, y, v*100, levels=cbar_levels, cmap='RdBu_r', extend='both')
csnoic = ax.contour( z, y, mach,  levels=[1.0], colors='lime', linewidths=2.0, zorder=9)

if streamline:
    ax.streamplot(z,y, w, v, color='black', linewidth=0.5, density=1.0)
    #ax.quiver(z[::4,::4], y[::4,::4], w[::4,::4], v[::4,::4], color='black', scale=100)

ax.fill_between( pz/length_unit, y_bottom, y2=ywall, color='gray', zorder=10 )
ax.set_aspect('equal')

cbar = plt.colorbar( cs, 
                     orientation='vertical', 
                     location='left', 
                     aspect=10,
                     ticks=cbar_ticks,
                     pad=0.30,
                     shrink=0.6) 

cbar.outline.set_linewidth(1.5)

cbar.ax.tick_params( direction='in',
                     left=True,right=False,
                     labelleft=True,labelright=False,
                     length=5,
                     width=1.0)

cbar.ax.set_xlabel(r'$\frac{\langle v \rangle}{u_{\infty}}\cdot 100$',labelpad=20)

ax.set_xlim( x_lim )
ax.set_ylim( y_lim )

ax.set_xticks( x_ticks )
ax.set_xticklabels( x_ticklabels )
ax.set_yticks( y_ticks )
ax.tick_params(which='major',
               axis='both', 
               direction='out',
               length=10.0,
               width=1.0)
ax.tick_params(axis='y', pad=15)
ax.tick_params(axis='x', pad=10)

ax.set_xlabel(x_label)
ax.set_ylabel(y_label)

ax.text( loc_tag[0], loc_tag[1],
         params.tag,
         va='center',
         ha='right',
         zorder=12,
         bbox={"fc":"white","alpha":0.8,"ec":"None"})    

ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(1.5)
ax.spines[:].set_zorder(11)


os.chdir( dirs.pp_statistics + '/yz_planes' )
plt.savefig(figname, dpi=300)
print(f"Figure saved to {dirs.pp_statistics}/yz_planes/{figname}.")