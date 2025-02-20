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
from   vista.params      import Params
from   vista.directories import Directories
from   vista.tools       import define_wall_shape
from   vista.tools       import lin_grow

# =============================================================================

From_File = True

# if From_file, give workpath to inca_grid.bin

workpath = '/home/wencan/temp/220927'
#workpath = '/media/wencanwu/Seagate Expansion Drive/temp/220825/results'
#workpath = '/home/wencanwu/my_simulation/STBLI_mid_Re'

# =============================================================================

os.chdir( workpath )

# plot from reading inca_grid.bin

if From_File:

    dirs           = Directories( workpath )
    gridfile       = dirs.grid
    parametersfile = dirs.case_para_file

    params   = Params( parametersfile )
    casecode = params.casecode
    delta    = params.delta_0

    # - read in grid info

    G = GridData( gridfile )
    G.read_grid()
    bbox = [-20.0,20.0,-5.0,40.0,-11.0,11.0]

    block_list, indx_slic = G.select_sliced_blockgrids( 'X', 0.1, bbox )

    fig, ax = plt.subplots( figsize = (20,24) )

    buff = 3

    for num in block_list:
        
        g = G.g[num-1]
        
        gz = (g.gz[buff:-buff+1] - 0.5*g.hz[buff:-buff+1])/delta
        gy = (g.gy[buff:-buff+1] - 0.5*g.hy[buff:-buff+1])/delta
        
        zmin = np.min( gz ) ; zmax = np.max( gz )
        ymin = np.min( gy ) ; ymax = np.max( gy )
        
        # draw horizontal lines
        for y in gy:
            ax.plot([zmin, zmax],[y,y],'black',linewidth=2)
        
        # draw vertical lines
        for z in gz:
            ax.plot([z, z],[ymin,ymax],'black',linewidth=2)
        
        if ymin < 0.0:
            y_wall = define_wall_shape(gz*delta,casecode=casecode,write=False)/delta
            
            ax.fill_between(gz,-0.2,y_wall,color='gray',zorder=10,alpha=0.7)

# plot from defining grid directly
else:
    
    delta = 5.2
    casecode = '241018'

    gy = np.arange(-1.254, 0.0, 0.01045)
    gy = np.concatenate( (gy,lin_grow(0.0, 0.01045, 1.0219,len=73)[0]) )
    
    gz = np.linspace(0,5.2,193)
    
    zmin = np.min( gz ) ; zmax = np.max( gz )
    ymin = np.min( gy ) ; ymax = np.max( gy )
    
    fig, ax = plt.subplots( figsize = (20,24) )
    
    # draw horizontal lines
    for y in gy:
        
        ax.plot([zmin/delta, zmax/delta],[y/delta,y/delta],'black',linewidth=1)
    
    # draw vertical lines
    for z in gz:
        
        ax.plot([z/delta, z/delta],[ymin/delta,ymax/delta],'black',linewidth=1)
    
    # block edges
    
    for y in [-1.254,-1.0032,-0.7524,-0.5016,-0.2508,0.0,0.325544991,0.873080888,1.7939845]:
        ax.plot([zmin/delta, zmax/delta],[y/delta,y/delta],'red',linewidth=1)
    
    for z in [0.0,1.3,2.6,3.9,5.2]:
        ax.plot([z/delta, z/delta],[ymin/delta,ymax/delta],'red',linewidth=1)
    
    # draw wall
    if ymin < 0.0:
        y_wall = define_wall_shape(gz,casecode=casecode,write=False)/delta
        
        ax.fill_between(gz/delta,-0.3,y_wall,color='gray',zorder=10,alpha=0.95)
    
ax.set_xlim([0,1])
ax.set_ylim([-0.2,1.0])
ax.set_aspect('equal', adjustable='box')
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

plt.show()

#plt.savefig("mesh.pdf")

plt.close()
    
    

    
    
    
    
