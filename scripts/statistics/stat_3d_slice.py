#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   stat_3d_slice.py
@Time    :   2026/01/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Slice 3D statistic data
'''

# need to install xvfbwrapper, and update 
# /path/to/conda/env/pp/lib/libstdc++.so.6 to have GLIBCXX_3.4.30

off_screen = False

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()
    
import os
import gc
import sys
import time
import numpy             as     np
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.mpi         import MPIenv
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.params      import Params
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.material    import get_visc
from   vista.directories import create_folder
from   vista.tools       import crop_border
from   vista.plot_setting import cpos_callback

def main():

    casefolder = '/home/wencan/temp/250710'
    bbox       = [-90.0, 60, -1.3, 31.0, 0.0, 999]
    clipbox    = [-25, 3, -1.0, 4.0, 0.0, 2.0]
    vars_read  = ['u','v','w','p','T']
    vars_out   = ['u','v','w','p','T','mach','grad_rho','DS']
    rescale    = [-50.4, 0.0, 0.0, 5.2, 5.2, 5.2]
    
    dirs       = Directories( casefolder )
    params     = Params( dirs.case_para_file )
    rough      = not ('smooth' in params.casecode.lower())
    
    dataset    = data_preparation( casefolder, bbox, vars_read, vars_out, rescale )
    visualize( dataset, rough, params, clipbox )

def visualize( dataset: pv.MultiBlock, rough:bool, params:Params, clipbox ):
    
    data = dataset.cell_data_to_point_data().combine()
    data = data.clip_box(clipbox, invert=False)
    
    if rough:
        wall = data.contour( [0.02], scalars='wd' )
    else:
        wall = data.slice( normal='Y', origin=(0,0.02,0) )
    
    data.set_active_scalars( 'DS' )
    bgslice = data.slice( normal='Z', origin=(0,0,0.01) )
    
    data['u'] = data['u'] / params.u_ref
    data.set_active_scalars( 'u' )
    xslices = []    
    x_locations = np.linspace( -23.0, -9, 8, endpoint=True)
    for x in x_locations:
        slicei = data.slice( normal='X', origin=(x,0,0) )
        xslices.append( slicei )
    
    bubble = data.contour( [-0.0001], scalars='u' )
    
    p = pv.Plotter( off_screen=off_screen, window_size=[1920,1080] )
    p.add_mesh( bgslice, cmap='gray', clim=[0.0, 0.5] )
    p.add_mesh( wall, color='gray' )
    for slicei in xslices:
        p.add_mesh( slicei, cmap='coolwarm', clim=[0.0, 1.0] )
    p.add_mesh( bubble, color='blue' )
    
    p.add_axes()
    p.set_background( 'white' )
    
    p.show()


def data_preparation( case_dir, bbox, vars_read, vars_out, rescale ):
    
    dirs     = Directories( case_dir )
    params   = Params( dirs.case_para_file )
    
    # - read in the grid data
    grid = GridData( dirs.grid )
    grid.read_grid()
    blocklist = grid.select_blockgrids( bbox, mode='within')

    if not 'smooth' in params.casecode:
        rough = True
        wdfile = get_filelist( dirs.wall_dist, 'snapshot.bin')[0]
        snapwd = Snapshot( wdfile )
        snapwd.read_snapshot( blocklist, ['wd'] )
    
    else:
        rough = False
    
    stat = StatisticData( dirs.statistics )
    stat.grid3d = grid
    stat.read_statistic( blocklist, vars_in=vars_read )
    stat.compute_vars(      blocklist, ['mach'] )
    stat.compute_gradients( blocklist, ['grad_rho','DS'] )
    
    # copy wall distance to statistic data
    if rough:
        for blnum in blocklist:
            stat.bl[stat.bl_nums.index(blnum)].df['wd'] = snapwd.snap_data[snapwd.bl_nums.index(blnum)].df['wd']
        vars_out += ['wd']
        
    dataset = pv.MultiBlock( stat.create_vtk_multiblock(blocklist, vars_out,buff=2, mode='oneside', rescale=rescale) )
    
    return dataset


# =============================================================================
if __name__ == "__main__":
    main()

