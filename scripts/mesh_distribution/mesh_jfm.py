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

# - define the method of drawing
def plot_block_edge( grid:GridData, block_list:list, ax, type, wall=False, lines=False ):

    for num in block_list:
        
        g  = grid.g[num-1]
        
        if type == 'xy':
            x0 = (g.lx0 - x_imp) / delta
            x1 = (g.lx1 - x_imp) / delta
        elif type == 'yz':
            x0 =  g.lz0          / delta
            x1 =  g.lz1          / delta
        
        y0 =  g.ly0          / delta
        y1 =  g.ly1          / delta
        
        if lines and type=='yz':
            
            gz = g.pz[3:-3] / delta
            gy = g.py[3:-3] / delta
            
            ax.vlines( gz[::4], y0, y1, color='gray', lw = 0.5)
            ax.hlines( gy[::4], x0, x1, color='gray', lw = 0.5)
        
        ax.plot( [x0,x1], [y0,y0], 'black', '-', lw=1 )
        ax.plot( [x1,x1], [y0,y1], 'black', '-', lw=1 )
        ax.plot( [x1,x0], [y1,y1], 'black', '-', lw=1 )
        ax.plot( [x0,x0], [y1,y0], 'black', '-', lw=1 )

    if wall:
        
        z_wall = np.linspace(-10.4,0.0,501)
        y_wall = define_wall_shape( z_wall, casecode=casecode, write=False ) / delta
        
        ax.fill_between(z_wall/delta, -0.3, y_wall, color='gray', zorder=10)
    
        
        

# =============================================================================
# - draw xy plane

ax  = fig.add_axes([0.05,0.20,0.70,0.4])
block_list, _ = grid.select_sliced_blockgrids( 'Z', 0.01 )

plot_block_edge( grid, block_list, ax, 'xy' )

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
ax.plot( [-13.78,-7.8], [1.95,8.6], 'blue', '-', lw=1 )

# -- large circle

sub_ax = fig.add_axes([0.35,0.4,0.18,0.18])

plot_block_edge( grid, block_list, sub_ax, 'xy' )

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

# =============================================================================
# - yz plane

ax2 = fig.add_axes([0.70,0.20,0.15,0.4])
block_list, _ = grid.select_sliced_blockgrids( 'X', -53.6 )

plot_block_edge( grid, block_list, ax2, 'yz' )

# -- set ranges
ax2.set_xlim([-3.0,3.0])
ax2.set_ylim([-1.0,17.0])

# -- hide top and right spines
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# -- set ticks
ax2.xaxis.set_major_locator( ticker.MultipleLocator(2)  )
ax2.set_yticks( np.arange(0,20,5) ) # some problem if using set_major_locator
ax2.set_yticklabels( [0,5,10,15]  ) # labels won't show
ax2.xaxis.set_minor_locator( ticker.MultipleLocator(1))
ax2.yaxis.set_minor_locator( ticker.MultipleLocator(2.5))
ax2.tick_params(which='major',
                axis='both', 
                direction='in',
                length=10.0,
                width=1.0, 
                labelsize=20,
                pad=10)
ax2.tick_params(which='minor',
                axis='both', 
                direction='in',
                length=5.0,
                width=1.0)
ax2.set_aspect('equal')

# -- set axis labels
ax2.set_xlabel(r'$z/\delta_0$')

# - draw the enlarged view as a inset subplot
# ----- small circle
circle2 = patches.Circle((-0.5,0.0),
                         radius=0.35,
                         fill=False,
                         edgecolor='blue')
ax2.add_patch( circle2 )

# ----- large circle

sub_ax2 = fig.add_axes([0.60,0.28,0.50,0.25])

plot_block_edge( grid, block_list, sub_ax2, 'yz', wall=True, lines=True )

circle3 = patches.Circle((0.5,0.5),
                         radius=0.5,
                         fill=True,
                         transform=sub_ax2.transAxes,
                         edgecolor='blue',
                         facecolor='white',
                         lw=2,
                         clip_on=True )

sub_ax2.add_patch( circle3 )
sub_ax2.set_clip_path(circle3)

for artist in sub_ax2.get_children():
    artist.set_clip_path(circle3)

sub_ax2.set_xlim([-0.7,-0.1])
sub_ax2.set_ylim([-0.20,0.40])
sub_ax2.set_aspect('equal')

sub_ax2.set_xticks([])
sub_ax2.set_yticks([])
sub_ax2.set_frame_on(False)

# ----- connection line

ax2.plot([-0.49,2.0],[0.22,6.0],'blue', lw=1)

# - save figure
# plt.savefig( '/home/wencan/temp/DataPost/midRe/mesh_config/mesh_jfm.pdf' )
plt.show()    