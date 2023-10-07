#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   draw_mesh.py
@Time    :   2023/10/07 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Script fro showing mesh and wall geometry
'''


import os

import sys

import numpy             as     np

import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData

from   vista.tools       import define_wall_shape
from   vista.tools       import read_case_parameter

# =============================================================================

workpath = '/home/wencanwu/my_simulation/temp/220927_lowRe/results'
workpath = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results'

# =============================================================================

os.chdir( workpath )

gridfile = workpath + '/inca_grid.bin'
parametersfile = workpath.split('/results')[0] + '/case_parameters'

parameters = read_case_parameter( parametersfile )
casecode = str( parameters.get('casecode') )
delta    = float( parameters.get('delta_0') )

# - read in grid info

G = GridData( gridfile )
G.read_grid()
bbox = [-10.0,10.0,-5.0,2.0, -11.0,11.0]

block_list, indx_slic = G.select_sliced_blockgrids( 'X', 0.0, bbox )

fig, ax = plt.subplots( figsize = (20,7) )

buff = 3

for num in block_list:
    
    g = G.g[num-1]
    
    gz = (g.gz[buff:-buff+1] - 0.5*g.hz[buff:-buff+1])/delta
    gy = (g.gy[buff:-buff+1] - 0.5*g.hy[buff:-buff+1])/delta
    
    
    zmin = np.min( gz ) ; zmax = np.max( gz )
    ymin = np.min( gy ) ; ymax = np.max( gy )
    
    # draw horizontal lines
    for y in gy:
        ax.plot([zmin, zmax],[y,y],'black',linewidth=0.8)
    
    # draw vertical lines
    for z in gz:
        ax.plot([z, z],[ymin,ymax],'black',linewidth=0.8 )
    
    y_wall = define_wall_shape(gz*delta, casecode=casecode, write=False)/delta
    
    ax.fill_between(gz,-0.2,y_wall,color='gray',zorder=10,alpha=0.7)
    
    
ax.set_xlim([0,1])
ax.set_ylim([-0.2,0.15])
ax.set_aspect('equal', adjustable='box')
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

#plt.show()

plt.savefig("mesh_"+casecode)

plt.close()
    
    

    
    
    
    
