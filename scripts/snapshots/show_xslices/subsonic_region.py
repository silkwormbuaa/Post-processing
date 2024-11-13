#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   subsonic_region.py
@Time    :   2024/11/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   shsow subsonic region change from all snapshots.
'''

off_screen = False

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()

import gc
import os
import sys
import time
import numpy              as     np
import pyvista            as     pv
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.mpi          import MPIenv
from   vista.grid         import GridData
from   vista.timer        import timer
from   vista.params       import Params
from   vista.snapshot     import Snapshot
from   vista.statistic    import StatisticData
from   vista.directories  import Directories
from   vista.tools        import crop_border
from   vista.tools        import get_filelist
from   vista.directories  import create_folder
from   vista.plot_setting import cpos_callback
from   vista.plot_setting import set_plt_rcparams

# - build MPI communication environment, set plotting parameters

mpi = MPIenv()

# =============================================================================

casedir  = '/home/wencan/temp/231124'

bbox      = [-30.0, 15.0, -1.3, 3.6, -999.0, 999.0]
vars_in   = ['u','v', 'w', 'T', 'p']
vars_out  = ['u','v','w','p','T','p_fluc','mach']

# =============================================================================

dirs       = Directories( casedir )
out_dir    = dirs.pp_snp_pfmax + '/figs'

params     = None
snapfiles  = None 
blocklist  = None
roughwall  = True
grid3d     = GridData()
stat       = None
snapwd     = Snapshot()

if mpi.rank == 0:
    
    create_folder( out_dir )
    
    params     = Params( dirs.case_para_file )
    roughwall  = params.roughwall

    snapfiles = get_filelist( dirs.snp_dir, 'snapshot.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot files.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    blocklist = grid3d.select_blockgrids( bbox, mode='within' )

    stat = StatisticData( dirs.statistics )
    stat.read_statistic( blocklist, ['p'] )

    if roughwall:
        snapwd = Snapshot( get_filelist(dirs.wall_dist, 'snapshot.bin')[0] )
        snapwd.read_snapshot( block_list=blocklist, var_read=['wd'] )
        print(f"Read in wall distance data from snapshot_{snapwd.itstep}.\n")

params     = mpi.comm.bcast( params,    root=0 )
roughwall  = mpi.comm.bcast( roughwall, root=0 )
snapfiles  = mpi.comm.bcast( snapfiles, root=0 )
blocklist  = mpi.comm.bcast( blocklist, root=0 )
grid3d     = mpi.comm.bcast( grid3d,    root=0 )
stat       = mpi.comm.bcast( stat,      root=0 )
snapwd     = mpi.comm.bcast( snapwd,    root=0 )
p_ref      = params.p_ref
u_ref      = params.u_ref

if roughwall: vars_out += ['wd']

x_imp      = params.x_imp
delta0     = params.delta_0
x_pfmax    = params.x_pfmax
x_pfmax    = x_pfmax*delta0 + x_imp

mpi.barrier()

os.chdir( out_dir )
clock = timer("show slice at pressure fluctuation max:")

def show_slice( snapfile ):

    snap        = Snapshot( snapfile )
    snap.grid3d = grid3d
    snap.read_snapshot( block_list=blocklist, var_read=vars_in )
    if roughwall:
        snap.copy_var_from( snapwd, ['wd'], blocklist=blocklist )
    
    snap.compute_vars( blocklist, ['mach'] )
    
    for bl_num in blocklist:
        
        snapblk = snap.snap_data[snap.bl_nums.index(bl_num)]
        statblk = stat.bl[stat.bl_nums.index(bl_num)]
        snapblk.df['p_fluc'] = snapblk.df['p'] - statblk.df['p']

    itstep  = snap.itstep
    itime   = snap.itime
    figname = f'slice_pfmax_{itstep:08d}.png'

    dataset = pv.MultiBlock(snap.create_vtk_multiblock( vars=vars_out, block_list=blocklist, mode='oneside'))
    sys.stdout.flush()
    
    dataset.set_active_scalars('u')
    point_data = dataset.cell_data_to_point_data().combine()

    uslicez     = point_data.slice(normal=[0,0,1], origin=[0,0,-10.3])

    if roughwall:
        point_data.set_active_scalars('wd')
        wallsurface = point_data.contour( [0.01] )
    else:
        wallsurface = point_data.slice( normal=[0.0,1.0,0.0], origin=[0.0,0.01,0.0] )
    
    wallsurface['p_fluc'] = wallsurface['p_fluc']/p_ref
    wallsurface.set_active_scalars('p_fluc')

    xslices = []
    x_slics = [x_pfmax-5.2,x_pfmax,x_pfmax+5.2]
    for j in range(len(x_slics)):
        xslices.append( point_data.slice( normal=[1,0,0], origin=[x_slics[j],0,0] ))    
    
    p = pv.Plotter(off_screen=off_screen,window_size=[1920,1080])
    p.add_mesh( uslicez, cmap='coolwarm', clim=[-0.8*u_ref,0.8*u_ref], 
                show_scalar_bar=False)
    
    p.add_mesh( wallsurface, cmap='coolwarm', 
                clim=[-0.4,0.4], show_scalar_bar=True)

    for slice in xslices:
        slice.set_active_scalars('u')
        vector2d        = np.zeros((slice.n_points,3))
        vector2d[:,1]   = slice['v']
        vector2d[:,2]   = slice['w']
        vector2d[::5,:] = 0.0
        slice['vector'] = vector2d
        sonline = slice.contour( [1.0], scalars='mach' )
        p.add_mesh( slice, cmap='coolwarm', clim=[-0.8*u_ref,0.8*u_ref], 
                    show_scalar_bar=False, lighting=False)
        p.add_arrows( slice.points, slice['vector'], mag=0.0015, color='black')
        p.add_mesh( sonline, color='green',  line_width=4.0 )
    
    camera_pos = [(x_pfmax-20.0,24.34,17.13),(x_pfmax,0.0,0.0),(0.52,0.80,-0.29)]
    p.camera_position = camera_pos
    cpos_callback( p )
    p.add_axes()
    p.show(figname)
    
    if off_screen:
        plt.figure(figsize=(6.4,3.6))
        plt.imshow(p.image)
        plt.plot([0.0,0.0],[10.0,10.0])
        plt.title(f"time={itime:6.2f} ms")
        plt.axis('off')                  # image axis, pixel count
        plt.tight_layout()
        plt.savefig(figname, dpi=300)
        plt.close()
    
    p.close()
    
    del snap, dataset, point_data, wallsurface, uslicez, xslices
    gc.collect()

# =============================================================================

if mpi.size == 1:
    print("No workers available. Master should do all tasks.")
    
    for i, snapfile in enumerate(snapfiles):
        show_slice( snapfile )
        clock.print_progress( i, len(snapfiles) )
        
else:
    if mpi.rank == 0:
        mpi.master_distribute( snapfiles )
    else:
        while True:
            task_index = mpi.worker_receive()
            
            if task_index is None: break
            else: 
                show_slice( snapfiles[task_index] )
                clock.print_progress( task_index, len(snapfiles), rank=mpi.rank )

mpi.barrier()

if mpi.rank == 0:

    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush() 
