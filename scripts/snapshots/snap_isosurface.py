#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   snap_isosurface.py
@Time    :   2024/08/26
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   visualise the isosurface of shock and vortices
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
import pyvista           as     pv
import matplotlib.pyplot as     plt
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.snapshot    import Snapshot
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.directories import create_folder

# - build MPI communication environment

comm    = MPI.COMM_WORLD
rank    = comm.Get_rank()
n_procs = comm.Get_size()

# =============================================================================
# option 
# =============================================================================

bbox      = [-30, 999, -1.0, 31.0, 5.2, 999]
gradients = ['Q_cr','div','vorticity','grad_rho','grad_rho_mod']
vars_out  = ['u','Q_cr','div','vorticity','grad_rho_mod']

snaps_dir = '/home/wencanwu/test/snapshots_231124'
gridfile  = '/media/wencan/Expansion/temp/231124/results/inca_grid.bin'
outdir    = '/home/wencanwu/test/snapshots_231124/output'

# =============================================================================

# - get snapshot file list/grid and distribute the tasks

snapfiles = None
grid3d    = GridData()

if rank == 0:
    
    snapfiles = get_filelist( snaps_dir, 'snapshot.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot files.")
    
    grid3d = GridData( gridfile )
    grid3d.read_grid()
    
snapfiles = comm.bcast( snapfiles, root=0 )
grid3d    = comm.bcast( grid3d,    root=0 )

n_snaps = len( snapfiles )
i_start, i_end = distribute_mpi_work(n_snaps, n_procs, rank)
snapfiles = snapfiles[i_start:i_end]

print(f"I am processor {rank:05d}, I take {len(snapfiles):5d} tasks.")
sys.stdout.flush()

comm.barrier()

# - read in snapshots and compute the gradients

os.chdir( create_folder(outdir) )
block_list = grid3d.select_blockgrids( bbox, mode='within' )
clock = timer("show isosurface")

for i,snapfile in enumerate(snapfiles):
    
    snap3d = Snapshot( snapfile )
    snap3d.verbose = False
    snap3d.grid3d = grid3d
    snap3d.read_snapshot( block_list )
    snap3d.compute_gradients( block_list, gradients )
#    snap3d.write_szplt( "test.szplt", vars=vars_out, block_list=block_list )

# -- generate the vtk dataset

    dataset = pv.MultiBlock(snap3d.create_vtk_multiblock( vars=vars_out, block_list=block_list, buff=3, mode='symmetry' ))
    sys.stdout.flush()

    dataset.set_active_scalars('u')
    uslicez = dataset.slice(normal=[0,0,1], origin=[0,0,5.3])
    uslicey = dataset.slice(normal=[0,1,0], origin=[0,0,0.05])
    
    point_data = dataset.cell_data_to_point_data().combine()

    sep_bubble = point_data.contour( [0.0] )
    
    point_data.set_active_scalars('grad_rho_mod')
    shock_front = point_data.contour( [0.2] )
    
    point_data.set_active_scalars('Q_cr')
    vortices = point_data.contour( [50000.0] )
    vortices.set_active_scalars('u')

# -- plot

# -- add figures elements
    
    # set the off_screen=True to avoid the interactive window
    # interactive window is not supported on remote servers, also there is a bug 
    # in vtk 9.x.x that the interactive window cannot be closed.
    
    p = pv.Plotter(off_screen=off_screen)
    cmapu = plt.get_cmap('RdBu_r',84)
    p.add_mesh(uslicez, opacity=1.0, clim=[-320,510],show_scalar_bar=True, cmap=cmapu)
    p.add_mesh(uslicey, opacity=1.0, clim=[-320,510],show_scalar_bar=True, cmap=cmapu)
    p.add_mesh(sep_bubble)
    p.add_mesh(shock_front, color='grey', opacity=0.5)
    p.add_mesh(vortices, cmap=cmapu, clim=[-320,510], show_scalar_bar=True)
    p.add_axes()    # add a vtk axes widget to the scene
    
# -- camera setting

    p.view_vector([0,0,0],viewup=[0.19,0.98,-0.18])
    p.camera.position = (-100,61,120)
    p.camera.focal_point = (37.8,8.8,-12.3)
    
# -- save the figure with matplotlib

    figname = f"isosurface_{snap3d.itstep:08d}"
    
    p.show(screenshot=figname + ".png")
    
    if off_screen:
        plt.imshow(p.image)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(figname + ".png", dpi=600)
        plt.close()

    p.close()
    
# - print the progress
    
    progress = (i+1)/len(snapfiles)
    print(f"Rank:{rank:05d},{i+1}/{len(snapfiles)} is done. " + clock.remainder(progress))
    print("------------------\n")
    
    del snap3d, dataset, p
    gc.collect()


if rank == 0:

    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush()       