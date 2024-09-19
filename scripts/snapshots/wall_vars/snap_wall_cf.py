#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   snap_wall_vars.py
@Time    :   2024/09/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   using pyvista to extract variables at the wall
'''

off_screen = False

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()
    
import os
import gc
import sys
import numpy             as     np
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.material    import get_visc

# =============================================================================

wd_file   = '/media/wencan/Expansion/temp/220927/supplements/wall_dist/snapshot_01329504/snapshot.bin'
snap_file = '/media/wencan/Expansion/temp/220927/supplements/wall_dist/snapshot_01329504/snapshot.bin'
grid_file = '/media/wencan/Expansion/temp/220927/results/inca_grid.bin'
bbox      = [-30.0,0.0,-3.0,2.0, 0.0,99.0]
vars_in   = ['u','T']

# =============================================================================

grid3d = GridData( grid_file )
grid3d.read_grid()
block_list = grid3d.select_blockgrids( bbox, mode='within' )

wd_snap = Snapshot( wd_file )
wd_snap.read_snapshot( block_list, var_read=['wd'] )

snap3d = Snapshot( snap_file )
snap3d.grid3d = grid3d
# snap3d.verbose = True
snap3d.read_snapshot( block_list=block_list, var_read=vars_in )

snap3d.copy_var_from( wd_snap, ['wd'] )

for bl in snap3d.snap_data:

    if bl.num in block_list:
        bl.df['mu'] = get_visc( np.array(bl.df['T']) )

# =============================================================================
# visualization
# =============================================================================

dataset = pv.MultiBlock( snap3d.create_vtk_multiblock(vars=vars_in,block_list=block_list,mode='oneside'))

point_data = dataset.cell_data_to_point_data().combine()
point_data.set_active_scalars('wd')

wallsurface = point_data.contour( [0.001] )
wallsurface.set_active_scalars('u')

print( wallsurface )

p = pv.Plotter(off_screen=off_screen)
cmap = plt.get_cmap('coolwarm',51)

p.add_mesh( wallsurface, cmap=cmap, show_scalar_bar=True)
p.add_axes()
p.show('test.png')

if off_screen:
    plt.imshow(p.image)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig("test.png", dpi=600)
    plt.close()

p.close()

gc.collect()


