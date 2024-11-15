#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   2d_zslice.py
@Time    :   2024/11/15 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   show shear layer vortices structure (relative velocity vector)
'''

off_screen = True

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()

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
from   vista.tools        import crop_to_contour_map
from   vista.tools        import get_filelist
from   vista.plane_analy  import compute_DS
from   vista.directories  import create_folder
from   vista.plot_setting import cpos_callback
from   vista.plot_setting import set_plt_rcparams

# - build MPI communication environment, set plotting parameters

mpi = MPIenv()
set_plt_rcparams(latex=False,fontsize=15)

# =============================================================================

casedir  = '/home/wencan/temp/231124'

vars_out = [ 'p_fluc' ]

varslist = ['u',                'T',               'p',               'DS',      
            'p_fluc',           'rho',            'rho_fluc',        'u_r',
            'v_r',              'w3',             'mach']
labels   = [r'$u/u_{\infty}$',     r'$T/T_{\infty}$',      r'$p/p_{\infty}$',         r'$DS$',   
            r"$p'/p_{\infty}$",    r"$\rho/\rho_{\infty}$",r"$\rho '/\rho_{\infty}$", r"$u_{r}/u_{\infty}$", 
            r"$v_{r}/u_{\infty}$", r'$\omega$',    r'$M$']
colmaps  = ['coolwarm',        'plasma',          'coolwarm',        'Greys_r', 
            'coolwarm',        'plasma',          'coolwarm',       'coolwarm',
            'coolwarm',        'coolwarm',        'coolwarm']
ranges   = [[-0.4,1.0],        [1.0,2.0],         [1.0,3.5],         [0.0,0.8], 
            [-0.5,0.5],        [0.5,2.5],         [-0.25,0.25],      [-0.2,0.2],        
            [-0.2,0.2],        [-1.5,1.5],        [0.0,2.0]]
vars_in  = ['u', 'v', 'w', 'T', 'p']
cutbox   = [-120.0, 120.0, -1.3, 86.0, 0.1, 0.11]
clipbox  = [-10, 5, 0, 4, -1, 1]

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
        create_folder( f'./{var}_vector/figs' )
    
    print( dirs.case_para_file)
    params = Params( dirs.case_para_file )
    
    snapfiles = get_filelist( dirs.snp_dir, 'snapshot_Z' )
    print(f"I am root, just found {len(snapfiles)} snapshot Z files.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    blocklist = grid3d.select_blockgrids( cutbox, mode='overlap' )
    
    statz = StatisticData( dirs.stat_zslice )
    statz.read_statistic(block_list=blocklist, vars_in=['u','v','p','rho'])
        
params    = mpi.comm.bcast( params,    root=0 )
snapfiles = mpi.comm.bcast( snapfiles, root=0 )
blocklist = mpi.comm.bcast( blocklist, root=0 )
grid3d    = mpi.comm.bcast( grid3d,    root=0 )
statz     = mpi.comm.bcast( statz,     root=0 )

u_ref     = params.u_ref
T_ref     = params.T_ref
p_ref     = params.p_ref
delta     = params.delta_0
rho_ref   = params.rho_ref
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
    snap.compute_gradients( grads=['grad_rho','vorticity'] )
    snap.compute_vars( blocklist, ['mach'] )

# -- data processing block by block

    for snapbl in snap.snap_data:
        
        bl_num = snapbl.num
        statbl = statz.bl[statz.bl_nums.index(bl_num)]
        
        snapbl.df['p_fluc']   = (snapbl.df['p']  - statbl.df['p']  )/p_ref
        snapbl.df['rho_fluc'] = (snapbl.df['rho']- statbl.df['rho'])/rho_ref
        snapbl.df['u_r']      = (snapbl.df['u']  - statbl.df['u']  )/u_ref
        snapbl.df['v_r']      = (snapbl.df['v']  - statbl.df['v']  )/u_ref
        snapbl.df['w3']       =  snapbl.df['w3'] /u_ref*delta
        snapbl.df['u']        =  snapbl.df['u']  /u_ref
        snapbl.df['v']        =  snapbl.df['v']  /u_ref
        snapbl.df['T']        =  snapbl.df['T']  /T_ref
        snapbl.df['p']        =  snapbl.df['p']  /p_ref
        snapbl.df['rho']      =  snapbl.df['rho']/rho_ref
        snapbl.df['DS']       = compute_DS( snapbl.df['grad_rho'], min=0.0, max=2.0)
        
# -- prepare for the visualization

    itstep      = snap.itstep
    itime       = snap.itime
    figname     = f"slice_z_{itstep:08d}.png"
    
    dataset     = pv.MultiBlock( snap.create_vtk_multiblock(rescale=rescale) )
    dataset     = dataset.cell_data_to_point_data().combine()
    dataset     = dataset.clip_box( clipbox, invert=False )
    sepline     = dataset.contour( [0.0], scalars='u' )
    sonline     = dataset.contour( [1.0], scalars='mach' )
    sys.stdout.flush()

# -- loop over all the output variables

    for var in vars_out:
        
        os.chdir( f'./{var}_vector/figs' )
        p = pv.Plotter( off_screen=off_screen, window_size=[1920,1080], border=False )

# ----- set color maps

        cmap    = plt.get_cmap(colmaps[varslist.index(var)],51)
        cmap.set_over('red')
        cmap.set_under('blue')
        
        dataset.set_active_scalars( var )
        p.add_mesh( dataset, 
                    cmap=cmap, 
                    clim=ranges[varslist.index(var)],
                    lighting=False, 
                    show_scalar_bar=False )
        
        p.add_mesh( sepline, color='yellow', line_width=4.0 )
        p.add_mesh( sonline, color='green',  line_width=4.0 )

# ----- add the relative velocity vector

        vector2d        = np.zeros((dataset.n_points,3))
        vector2d[:,0]   = dataset['u_r']
        vector2d[:,1]   = dataset['v_r']
        mask            = np.zeros_like( vector2d, dtype=bool )
        mask[::7,::]    = True
        vector2d        = vector2d*mask
        dataset['vector'] = vector2d
        p.add_arrows( dataset.points, dataset['vector'], mag=0.7, color='black')

# ----- set the view and render the image

        p.view_vector([0.0,0.0,1.0],viewup=[0.0,1.0,0.0])
        p.camera.tight()
        image = p.screenshot(return_img=True)
        p.close()

# ----- pass the image from pyvista to matplotlib

#        image = crop_to_contour_map(image)

        fig, ax = plt.subplots(figsize=(12.8,7.2))

        img = ax.imshow(image, extent=clipbox[:4], cmap=cmap, clim=ranges[varslist.index(var)])

        ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
        ax.set_ylabel(r'$y/\delta_0$')

        cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5 ) 
        cbar.ax.set_ylabel( labels[varslist.index(var)], loc='center')

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