#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   mesh_jfm.py
@Time    :   2025/02/20 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Illustration of block and mesh distributions.
'''


import os
import sys
import numpy              as     np
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker
import matplotlib.patches as     patches

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid         import GridData
from   vista.params       import Params
from   vista.directories  import Directories
from   vista.tools        import define_wall_shape

plt.rcParams["text.usetex"]         = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family']         = "Times New Roman"
plt.rcParams['font.size']           = 30

# =============================================================================

casefolder = '/home/wencan/temp/231124'

# =============================================================================
# - Read grid data
dirs       = Directories( casefolder )
grid       = GridData( dirs.grid )
params     = Params( dirs.case_para_file )
casecode   = params.casecode
delta      = params.delta_0
x_imp      = params.x_imp

# - read in grid info

grid.read_grid()

fig = plt.figure( figsize=(15,12)  )
ax  = fig.add_axes([0.10,0.20,0.9,0.6])

# - define the method of drawing
def plot_block_edge( grid:GridData, block_list:list, ax ):

    for num in block_list:
        
        g  = grid.g[num-1]
        x0 = (g.lx0 - x_imp) / delta
        x1 = (g.lx1 - x_imp) / delta
        y0 =  g.ly0          / delta
        y1 =  g.ly1          / delta
        
        ax.plot( [x0,x1], [y0,y0], 'black', '-', lw=1 )
        ax.plot( [x1,x1], [y0,y1], 'black', '-', lw=1 )
        ax.plot( [x1,x0], [y1,y1], 'black', '-', lw=1 )
        ax.plot( [x0,x0], [y1,y0], 'black', '-', lw=1 )
        

# - draw xy plane

block_list, _ = grid.select_sliced_blockgrids( 'Z', 0.01 )

plot_block_edge( grid, block_list, ax )

# -- set ranges
ax.set_xlim([-33.0,15.0])
ax.set_ylim([-1.0 ,17.0])

# -- hide top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# -- set ticks
ax.xaxis.set_major_locator( ticker.MultipleLocator(5)  )
ax.set_yticks( np.arange(0,20,5) ) # some problem if using set_major_locator
ax.set_yticklabels( [0,5,10,15]  ) # labels won't show
ax.xaxis.set_minor_locator( ticker.MultipleLocator(2.5))
ax.yaxis.set_minor_locator( ticker.MultipleLocator(2.5))
ax.tick_params(which='major',
               axis='both', 
               direction='in',
               length=10.0,
               width=1.0, 
               labelsize=20,
               pad=10)
ax.tick_params(which='minor',
               axis='both', 
               direction='in',
               length=5.0,
               width=1.0)
ax.set_aspect('equal')

# -- set axis labels
ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
ax.set_ylabel(r'$y/\delta_0$')

# - draw the enlarged view as a inset subplot
# -- small circle
circle = patches.Circle((-15,1.0),
                        radius=1.5,
                        fill=False,
                        edgecolor='blue')
ax.add_patch( circle )
ax.plot( [-13.78,-7.0], [1.95,8.0], 'blue', '-', lw=1 )

# -- large circle

sub_ax = fig.add_axes([0.5,0.45,0.24,0.30])

plot_block_edge( grid, block_list, sub_ax )

circle1 = patches.Circle((0.5,0.5), 
                         radius=0.5, 
                         fill=True,
                         transform=sub_ax.transAxes,
                         edgecolor='blue', 
                         facecolor='white',
                         lw=2,
                         clip_on=True )

sub_ax.add_patch( circle1 )
sub_ax.set_clip_path(circle1)

for artist in sub_ax.get_children():
    artist.set_clip_path(circle1)

sub_ax.set_xlim([ -16.5,-13.5])
sub_ax.set_ylim([-0.5,2.5])
sub_ax.set_aspect('equal')

# -- hide frames of sub_ax
sub_ax.set_xticks([])
sub_ax.set_yticks([])
sub_ax.set_frame_on(False)

plt.show()    