#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   snap_wall_vars.py
@Time    :   2024/09/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   using pyvista to visualise cf at the wall
'''

# need to install xvfbwrapper, and update 
# /path/to/conda/env/pp/lib/libstdc++.so.6 to have GLIBCXX_3.4.30

off_screen = True

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()
    
import os
import gc
import sys
import time
import numpy              as     np
import pyvista            as     pv
import matplotlib.pyplot  as     plt
from   mpi4py             import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid         import GridData
from   vista.timer        import timer
from   vista.snapshot     import Snapshot
from   vista.directories  import Directories
from   vista.tools        import get_filelist
from   vista.tools        import read_case_parameter
from   vista.tools        import distribute_mpi_work
from   vista.material     import get_visc
from   vista.directories  import create_folder
from   vista.plot_setting import cpos_callback

# - build MPI communication environment

comm    = MPI.COMM_WORLD
rank    = comm.Get_rank()
n_procs = comm.Get_size()

# =============================================================================
# option
# =============================================================================

case_dir  = '/home/wencanwu/test/220927'
bbox      = [-30.0,110.0,-3.0, 0.01, -99.0,99.0]
walldist  = 0.019
vars_in   = ['u','T']

# =============================================================================
# preparation
# =============================================================================

dirs       = Directories( case_dir )
out_dir    = dirs.pp_snp_fricpv

p_dyn      = None
snapfiles  = None
block_list = None
roughwall  = True
params     = None
grid3d     = GridData()
wd_snap    = Snapshot()

if rank == 0:
    
    create_folder( out_dir )
    
    params    = read_case_parameter( dirs.case_para_file )
    u_ref     = float(params.get('u_ref'))
    rho_ref   = float(params.get('rho_ref'))
    p_dyn     = 0.5 * rho_ref * u_ref**2
    roughwall = True if str(params.get('roughwall')).lower() == 'true' else False
    
    snapfiles = get_filelist( dirs.snp_dir, 'snapshot.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot files.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    block_list = grid3d.select_blockgrids( bbox, mode='within' )
    
    if roughwall:
        wd_file = get_filelist( dirs.wall_dist, 'snapshot.bin' )[0]
        wd_snap = Snapshot( wd_file )
        wd_snap.read_snapshot( block_list, var_read=['wd'] )
        print(f"Wall distance from {wd_file} has been read.\n")

params     = comm.bcast( params,     root=0 )
p_dyn      = comm.bcast( p_dyn,      root=0 )
roughwall  = comm.bcast( roughwall,  root=0 )
snapfiles  = comm.bcast( snapfiles,  root=0 )
block_list = comm.bcast( block_list, root=0 )
grid3d     = comm.bcast( grid3d,     root=0 )
wd_snap    = comm.bcast( wd_snap,    root=0 )
Re_ref     = float(params.get('Re_ref'))
visc_law   = params.get('visc_law')

vars_out   = vars_in
if roughwall : vars_out += ['mu','wd']
else         : vars_out += ['mu']

n_snaps   = len( snapfiles )
i_s, i_e  = distribute_mpi_work(n_snaps, n_procs, rank)
snapfiles = snapfiles[i_s:i_e]

print(f"I am processor {rank:05d}, I take {len(snapfiles):5d} tasks.")
sys.stdout.flush()

comm.barrier()

os.chdir( out_dir )
clock = timer("show cf")

# -- loop over the snapshots

for i, snap_file in enumerate(snapfiles):

    snap3d = Snapshot( snap_file )
    snap3d.grid3d = grid3d
    snap3d.read_snapshot( block_list=block_list, var_read=vars_in )

    itstep  = snap3d.itstep
    itime   = snap3d.itime
    figname = f'cf_{itstep:08d}.png'
    
    if roughwall:
        snap3d.copy_var_from( wd_snap, ['wd'], block_list )

    for bl in snap3d.snap_data:

        if bl.num in block_list:
            bl.df['mu'] = get_visc( np.array(bl.df['T']), Re_ref, law=visc_law )

# =============================================================================
# visualization
# =============================================================================
    
    dataset    = pv.MultiBlock( snap3d.create_vtk_multiblock(vars=vars_out,block_list=block_list,mode='oneside') )
    point_data = dataset.cell_data_to_point_data().combine()
    
    if roughwall:
        point_data.set_active_scalars('wd')
        wallsurface = point_data.contour( [walldist] )
        friction    = wallsurface['mu']*wallsurface['u']/wallsurface['wd']

    else:
        wallsurface = point_data.slice( normal=[0.0,1.0,0.0], origin=[0.0,walldist,0.0] )
        friction    = wallsurface['mu']*wallsurface['u']/walldist

    wallsurface['cf'] = friction / p_dyn
    wallsurface.set_active_scalars('cf')

    p = pv.Plotter(off_screen=off_screen, window_size=[1920,1080])
    cmap = plt.get_cmap('RdBu_r',51)

    p.add_mesh( wallsurface,    cmap=cmap, clim=[-0.006,0.006], show_scalar_bar=True, 
                lighting=False)
    
    p.view_vector([0.0,1.0,0.0],viewup=[0.0,0.0,-1.0])
    camera_pos = [(40.0,150.0,0.0),(40.0,0.0,0.0),(0.0,0.0,-1.0)]
    p.camera_position = camera_pos
#    cpos_callback( p )
    
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

    del snap3d,p,dataset,point_data,wallsurface,friction
    gc.collect()

# - print the progress
    
    progress = (i+1)/len(snapfiles)
    print(f"Rank:{rank:05d},{i+1}/{len(snapfiles)} is done. " + clock.remainder(progress))
    print("------------------\n")
    sys.stdout.flush()

comm.barrier()

if rank == 0:

    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush() 
