#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   2d_zslice.py
@Time    :   2024/10/24 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   read in snapshot_Z_xxxx.bin and show it with pyvista and matplotlib
'''

off_screen = True

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()

import os
import sys
import time
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
from   vista.plane_analy  import compute_DS
from   vista.directories  import create_folder
from   vista.plot_setting import cpos_callback
from   vista.plot_setting import set_plt_rcparams

# - build MPI communication environment, set plotting parameters

mpi = MPIenv()
set_plt_rcparams(latex=False,fontsize=15)

# =============================================================================

casedir  = '/home/wencan/temp/smooth_mid'
vars_in  = ['u', 'T', 'p']
vars_out = ['u', 'T', 'p', 'DS','p_fluc']
labels   = [r'$u/u_{\infty}$', r'$T/T_{\infty}$', r'$p/p_{\infty}$', r'$DS$',   r"$p'/p_{\infty}$"]
colmaps  = ['coolwarm',        'plasma',          'coolwarm',        'Greys_r', 'coolwarm']
ranges   = [[-0.4,1.0],        [1.0,2.0],         [1.0,3.5],         [0.0,0.8], [-0.5,0.5]]
cutbox   = [-120.0, 120.0, -1.3, 86.0, 0.1, 0.11]
clipbox  = [-20, 12, 0, 10, -1, 1]

# =============================================================================

dirs     = Directories( casedir )
out_dir  = dirs.pp_snp_zslice

# allocate variable that will be broadcasted to all workers

params    = None
snapfiles = None
blocklist = None
grid3d    = None
statz     = None

# root does the preparation

if mpi.is_root:
    
    os.chdir( create_folder( out_dir ) )
    
    for var in vars_out:
        create_folder( f'./{var}/figs' )
    
    print( dirs.case_para_file)
    params = Params( dirs.case_para_file )
    
    snapfiles = get_filelist( dirs.snp_dir, 'snapshot_Z' )
    print(f"I am root, just found {len(snapfiles)} snapshot Z files.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    blocklist = grid3d.select_blockgrids( cutbox, mode='overlap' )
    
    statz = StatisticData( dirs.stat_zslice )
    statz.read_statistic(block_list=blocklist, vars_in=['p'])
        
params    = mpi.comm.bcast( params,    root=0 )
snapfiles = mpi.comm.bcast( snapfiles, root=0 )
blocklist = mpi.comm.bcast( blocklist, root=0 )
grid3d    = mpi.comm.bcast( grid3d,    root=0 )
statz     = mpi.comm.bcast( statz,     root=0 )

u_ref     = params.u_ref
T_ref     = params.T_ref
p_ref     = params.p_ref
delta     = params.delta_0
rescale   = [-params.x_imp, 0.0, 0.0, delta, delta, delta]

mpi.barrier()

os.chdir( out_dir )
clock = timer("show slices from snapshot_Z_xxxx.bin:")

# read in the snapshots and show them

def show_slice( snapfile ):

# -- read in the snapshot and compute density gradient

    snap        = Snapshot( snapfile )
    snap.grid3d = grid3d
    snap.read_snapshot( block_list=blocklist, var_read=vars_in )
    snap.compute_gradients( grads=['grad_rho'] )

# -- data processing block by block

    for snapbl in snap.snap_data:
        
        bl_num = snapbl.num
        statbl = statz.bl[statz.bl_nums.index(bl_num)]
        snapbl.df['p_fluc'] = (snapbl.df['p']-statbl.df['p'])/p_ref
        
        snapbl.df['u'] = snapbl.df['u']/u_ref
        snapbl.df['T'] = snapbl.df['T']/T_ref
        snapbl.df['p'] = snapbl.df['p']/p_ref
        snapbl.df['DS']= compute_DS( snapbl.df['grad_rho'], min=0.0, max=2.0)
        
# -- prepare for the visualization

    itstep      = snap.itstep
    itime       = snap.itime
    figname     = f"slice_z_{itstep:08d}.png"
    
    dataset     = pv.MultiBlock( snap.create_vtk_multiblock(rescale=rescale) )
    dataset     = dataset.cell_data_to_point_data().combine()
    dataset     = dataset.clip_box( clipbox, invert=False )
    sepline     = dataset.contour( [0.0], scalars='u' )
    sys.stdout.flush()

# -- loop over all the output variables

    for i, var in enumerate(vars_out):
        
        os.chdir( f'./{var}/figs' )
        p = pv.Plotter( off_screen=off_screen, window_size=[1920,1080], border=False )

# ----- set color maps

        cmap    = plt.get_cmap(colmaps[i],51)
        cmap.set_over('red')
        cmap.set_under('blue')
        
        dataset.set_active_scalars( var )
        p.add_mesh( dataset, 
                    cmap=cmap, 
                    clim=ranges[i],
                    lighting=False, 
                    show_scalar_bar=False )
        
        p.add_mesh( sepline, color='yellow', line_width=4.0 )

        p.view_vector([0.0,0.0,1.0],viewup=[0.0,1.0,0.0])
        p.camera.tight()
        image = p.screenshot(return_img=True)
        p.close()

# ----- pass the image from pyvista to matplotlib

        image = crop_border(image)

        fig, ax = plt.subplots(figsize=(12.8,7.2))

        img = ax.imshow(image, extent=clipbox[:4], cmap=cmap, clim=ranges[i])

        ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
        ax.set_ylabel(r'$y/\delta_0$')

        cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5 ) 
        cbar.ax.set_ylabel( labels[i], loc='center')

        plt.title(f"t = {itime:.2f} ms", loc='center',y=1.05)

        if off_screen:
            plt.savefig( figname, dpi=150)
        else:
            plt.show()
        
        plt.close()
        os.chdir( '../../' )

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