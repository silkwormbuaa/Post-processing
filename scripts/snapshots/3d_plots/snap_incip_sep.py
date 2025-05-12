#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   snap_pressure_fluc.py
@Time    :   2025/04/28 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Show the reversed flow isosurface at the wall
'''

# need to install xvfbwrapper, and update 
# /path/to/conda/env/pp/lib/libstdc++.so.6 to have GLIBCXX_3.4.30

off_screen = True

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
from   vista.directories  import Directories
from   vista.tools        import get_filelist
from   vista.material     import get_visc
from   vista.directories  import create_folder
from   vista.plot_setting import cpos_callback

# - build MPI communication environment

mpi = MPIenv()

# =============================================================================
# option
# =============================================================================

casedir   = '/home/wencan/temp/241030/'
bbox      = [-35.0, 40.0, -1.3, 20, -999.0, 999.0]
vars_in   = ['u','p','T']
vars_out  = ['u','p','T']
#gradients = ['Q_cr','div','vorticity','grad_rho','grad_rho_mod']


# =============================================================================

dirs       = Directories( casedir )
out_dir    = dirs.pp_snp_incipsep + '/figs'

params     = None
p_dyn      = None
snapfiles  = None 
blocklist  = None
roughwall  = True
grid3d     = GridData()
snapwd     = Snapshot()

if mpi.rank == 0:
    
    create_folder( out_dir )
    
    params     = Params( dirs.case_para_file )
    roughwall  = params.roughwall
    u_ref     = params.u_ref
    rho_ref   = params.rho_ref
    p_dyn     = 0.5 * rho_ref * u_ref**2

    snapfiles = get_filelist( dirs.snp_dir, 'snapshot.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot files.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    blocklist = grid3d.select_blockgrids( bbox, mode='within' )

    if roughwall:
        snapwd = Snapshot( get_filelist(dirs.wall_dist, 'snapshot.bin')[0] )
        snapwd.read_snapshot( block_list=blocklist, var_read=['wd'] )
        print(f"Read in wall distance data from snapshot_{snapwd.itstep}.\n")

params     = mpi.comm.bcast( params,    root=0 )
p_dyn      = mpi.comm.bcast( p_dyn,     root=0 )
roughwall  = mpi.comm.bcast( roughwall, root=0 )
snapfiles  = mpi.comm.bcast( snapfiles, root=0 )
blocklist  = mpi.comm.bcast( blocklist, root=0 )
grid3d     = mpi.comm.bcast( grid3d,    root=0 )
snapwd     = mpi.comm.bcast( snapwd,    root=0 )
p_ref      = params.p_ref
u_ref      = params.u_ref
Re_ref     = params.Re_ref
visc_law   = params.visc_law

casecode   = params.casecode
if casecode in ['smooth_mid','241030','231124','241018']:
    walldist = 0.005
else:
    walldist = 0.019

if roughwall: vars_out += ['mu','wd']
else:         vars_out += ['mu']

x_imp      = params.x_imp
delta0     = params.delta_0
x_incip    = params.x_incip
x_incip    = x_incip*delta0 + x_imp

mpi.barrier()

os.chdir( out_dir )
clock = timer("show reversed flow isosurface near pressure fluctuation max:")

def show_slice( snapfile ):

    snap        = Snapshot( snapfile )
    snap.grid3d = grid3d
    snap.read_snapshot( block_list=blocklist, var_read=vars_in )
    if roughwall:
        snap.copy_var_from( snapwd, ['wd'], blocklist=blocklist )

    for bl in snap.snap_data:
        if bl.num in blocklist:
            bl.df['mu'] = get_visc( np.array(bl.df['T']), Re_ref, law=visc_law )
            
    itstep  = snap.itstep
    itime   = snap.itime
    figname = f'incip_sep_{itstep:08d}.png'

    dataset = pv.MultiBlock(snap.create_vtk_multiblock( vars=vars_out, block_list=blocklist, mode='oneside'))
    sys.stdout.flush()
    
    dataset.set_active_scalars('p')
    point_data = dataset.cell_data_to_point_data().combine()

    pslicez    = point_data.slice(normal=[0,0,1], origin=[0,0,-10.3])

    if roughwall:
        point_data.set_active_scalars('wd')
        wallsurface = point_data.contour( [walldist] )
        friction    = wallsurface['mu']*wallsurface['u']/wallsurface['wd']
    else:
        wallsurface = point_data.slice( normal=[0.0,1.0,0.0], origin=[0.0,walldist,0.0] )
        friction    = wallsurface['mu']*wallsurface['u']/walldist

    wallsurface['cf'] = friction/p_dyn
    wallsurface.set_active_scalars('cf')
    
    ruler_ticks = np.linspace(-15,-5,41, endpoint=True)
    points = list()
    for x_norm in ruler_ticks:
        x = x_norm * delta0 + x_imp
        if x > bbox[0] and x < bbox[1]:
            point1 = [x, walldist, -2.0*delta0]
            point2 = [x, walldist,  2.0*delta0]
            points.append( point1 )
            points.append( point2 )
    
    points = np.array(points)
    
    bubble = point_data.contour( [0.0], scalars='u' )
    bubble.set_active_scalars('p')
    
    p = pv.Plotter(off_screen=off_screen,window_size=[1920,1080])
    p.add_mesh( pslicez, cmap='coolwarm', clim=[40000,90000], 
                show_scalar_bar=False)

    cmap = plt.get_cmap('RdBu_r',51)
    p.add_mesh( wallsurface, cmap=cmap, 
                clim=[-0.007,0.007], show_scalar_bar=True,lighting=False)
    p.add_lines( points, color='black', width=1.0 )

    cmap = plt.get_cmap('coolwarm',51)
    p.add_mesh( bubble, cmap=cmap, clim=[40000,90000], show_scalar_bar=True,
                lighting=True)
    
    camera_pos = [(x_incip-25.0,31.0,25.0),(x_incip+10,0.0,0.0),(0.52,0.79,-0.33)]
    p.camera_position = camera_pos
    # cpos_callback( p )
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
    
    del snap, dataset, point_data, wallsurface
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