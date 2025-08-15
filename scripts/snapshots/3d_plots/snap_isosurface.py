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

off_screen = True

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
from   vista.params      import Params
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.material    import get_visc
from   vista.directories import create_folder
from   vista.tools       import crop_border
from   vista.plot_setting import cpos_callback

# - build MPI communication environment

mpi = MPIenv()

# =============================================================================
# option 
# =============================================================================

casefolder = '/home/wencan/temp/250710'

bbox      = [-58.0, 999, -1.3, 31.0, -999, 999]
gradients = ['Q_cr','div','vorticity','grad_rho','grad_rho_mod']
vars_out  = ['u','Q_cr','grad_rho_mod','p']

dirs      = Directories( casefolder )
snaps_dir = dirs.snp_dir
gridfile  = dirs.grid
outdir    = dirs.pp_snp_3dshock + '/figs'

# =============================================================================

# - get snapshot file list/grid and distribute the tasks

params     = None
p_dyn      = None
snapfiles  = None
block_list = None
snapwd     = None
grid3d     = GridData()

if mpi.is_root:
    
    snapfiles = get_filelist( snaps_dir, 'snapshot.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot files.")
    
    grid3d = GridData( gridfile )
    grid3d.read_grid()
    
    block_list = grid3d.select_blockgrids( bbox, mode='within' )
    
    params = Params( dirs.case_para_file )
    u_ref     = params.u_ref
    rho_ref   = params.rho_ref
    p_dyn     = 0.5 * rho_ref * u_ref**2
    
    if params.roughwall:
    
        wdfile = get_filelist( dirs.wall_dist, 'snapshot.bin' )[0]
        snapwd = Snapshot( wdfile )
        snapwd.read_snapshot(block_list, ['wd'])
    
    create_folder(outdir)
    
snapfiles  = mpi.comm.bcast( snapfiles,  root=0 )
grid3d     = mpi.comm.bcast( grid3d,     root=0 )
block_list = mpi.comm.bcast( block_list, root=0 )
snapwd     = mpi.comm.bcast( snapwd,     root=0 )
params     = mpi.comm.bcast( params,     root=0 )
p_dyn      = mpi.comm.bcast( p_dyn,      root=0 )
roughwall  = params.roughwall
Re_ref     = params.Re_ref
visc_law   = params.visc_law
highRe     = True if Re_ref > 9000 else False

if highRe: 
    walldist = 0.005
    grad_rho_limit = 0.2
else: 
    walldist = 0.019     
    grad_rho_limit = 0.15

if roughwall: vars_out += ['mu','wd']
else:         vars_out += ['mu']

sys.stdout.flush()
mpi.comm.barrier()

# - read in snapshots and compute the gradients

os.chdir( outdir )
clock = timer("show isosurface")

def plot_isosurface( snapfile ):
    
    snap3d = Snapshot( snapfile )
    snap3d.verbose = False
    snap3d.grid3d = grid3d
    snap3d.read_snapshot( block_list )
    if roughwall:
        snap3d.copy_var_from( snapwd, ['wd'], block_list )

    for bl in snap3d.snap_data:
        if bl.num in block_list:
            bl.df['mu'] = get_visc( np.array(bl.df['T']), Re_ref, law=visc_law )

    snap3d.compute_gradients( block_list, gradients )
#    snap3d.write_szplt( "test.szplt", vars=vars_out, block_list=block_list )

# -- generate the vtk dataset

    dataset = pv.MultiBlock(snap3d.create_vtk_multiblock( vars=vars_out, block_list=block_list, buff=3, mode='symmetry' ))
    
    # delete the unneeded variables to save memory
    itstep = snap3d.itstep
    itime  = snap3d.itime
    del snap3d ; gc.collect()
    sys.stdout.flush()
    
    point_data = dataset.cell_data_to_point_data().combine()
    point_data['p'] = point_data['p']/params.p_ref
    point_data['u'] = point_data['u']/params.u_ref
    point_data.set_active_scalars('p')
    pslicez = point_data.slice(normal=[0,0,1], origin=[0,0,-10.3])
    
    sep_bubble = point_data.contour(  [0.0], scalars='u' )
    shock_front = point_data.contour( [grad_rho_limit], scalars='grad_rho_mod' )
    
    vortices = point_data.contour( [50000.0], scalars='Q_cr' )
    vortices.set_active_scalars('u')

    if roughwall:
        wallsurface = point_data.contour( [walldist], scalars='wd' )
    else:
        wallsurface = point_data.slice( normal=[0.0,1.0,0.0], origin=[0.0,walldist,0.0] )

#    friction    = wallsurface['mu']*wallsurface['u']*params.u_ref/walldist
#    wallsurface['cf'] = friction/p_dyn
#    wallsurface.set_active_scalars('cf')

# -- plot

# -- add figures elements
    
    # set the off_screen=True to avoid the interactive window
    # interactive window is not supported on remote servers, also there is a bug 
    # in vtk 9.x.x that the interactive window cannot be closed.
    
    p = pv.Plotter(off_screen=off_screen,window_size=[1920,1080])
    
    cmapp = plt.get_cmap('coolwarm',51)
    p.add_mesh(pslicez,     opacity=1.0, clim=[1.0,3.5], show_scalar_bar=True, cmap=cmapp)
    
    cmapcf = plt.get_cmap('RdBu_r',51)
    p.add_mesh(wallsurface, opacity=1.0, color='gray', lighting=True )
    
    p.add_mesh(sep_bubble)
    p.add_mesh(shock_front, color='grey', opacity=0.5)
    p.add_mesh(vortices, cmap=cmapp, clim=[-0.2,1.0], show_scalar_bar=True)
    p.add_axes()    # add a vtk axes widget to the scene
    
# -- camera setting

    p.view_vector([0,0,0],viewup=[0.39,0.90,-0.20])
    p.camera.position = (-120,50,25)
    p.camera.focal_point = (20,-18,-5)
    
    # cpos_callback( p )
    # p.show()
    
# -- save the figure with matplotlib

    figname = f"isosurface_{itstep:08d}"
    
    image = p.screenshot(return_img=True)
    image = crop_border(image)
    
    if off_screen:
        fig = plt.figure(figsize=(6.4,3.6))
        # set the scope of the image in plt
        fig.subplots_adjust(left=0.0, right=1.0, top=0.9, bottom=0.0)
        
        ax  = fig.add_subplot(111)
        img = ax.imshow(image)
        ax.axis('off')
        
        plt.title(f"time = {itime:.2f} ms")
        plt.tight_layout()
        plt.savefig(figname + ".png", dpi=300)
        plt.close()

    p.close()
    
# - print the progress
    
    del dataset, p
    gc.collect()


if mpi.size == 1:
    print("No workers available. Master should do all tasks.")
    
    for i, snapfile in enumerate(snapfiles):
        plot_isosurface( snapfile )
        clock.print_progress( i, len(snapfiles) )
        
else:
    if mpi.rank == 0:
        mpi.master_distribute( snapfiles )
    else:
        while True:
            task_index = mpi.worker_receive()
            
            if task_index is None: break
            else: 
                plot_isosurface( snapfiles[task_index] )
                clock.print_progress( task_index, len(snapfiles), rank=mpi.rank )

mpi.barrier()

if mpi.rank == 0:

    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush() 