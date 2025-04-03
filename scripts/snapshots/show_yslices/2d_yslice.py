#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   2d_yslice.py
@Time    :   2025/02/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Read in snapshot_Y_xxxx.bin files and show it with pyvista and matplotlib
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
from   vista.tools        import get_filelist
from   vista.tools        import define_wall_shape
from   vista.plane_analy  import compute_DS
from   vista.directories  import create_folder
from   vista.plot_setting import cpos_callback
from   vista.plot_setting import set_plt_rcparams

# - build MPI communication environment, set plotting parameters

mpi = MPIenv()
set_plt_rcparams(latex=False,fontsize=15)

# =============================================================================

casedir  = '/home/wencan/temp/smooth_mid/'

locs     = [0.025]
locs     = np.array(locs) * 5.2

vars_out = ['u',                'v',               'w',                'p',
            'T',                'rho',             'mach',
            'u_r',              'v_r',             'w_r',              'p_fluc',
            'rho_fluc']

varslist = ['u',                'v',               'w',                'p',
            'T',                'rho',             'mach',
            'u_r',              'v_r',             'w_r',              'p_fluc',
            'rho_fluc']

labels   = [r'$u/u_{\infty}$',        r'$v/u_{\infty}$',       r'$w/u_{\infty}$',     r'$p/p_{\infty}$',
            r'$T/T_{\infty}$',        r'$\rho/\rho_{\infty}$', r'$Mach$',
            r'$u_{r}/u_{\infty}$',    r'$v_{r}/u_{\infty}$',   r'$w_{r}/u_{\infty}$', r"$p'/p_{\infty}$",
            r"$\rho'/\rho_{\infty}$"]

colmaps  = ['coolwarm',        'coolwarm',        'coolwarm',        'coolwarm',
            'plasma',          'plasma',          'coolwarm',
            'coolwarm',        'coolwarm',        'coolwarm',        'coolwarm',
            'coolwarm']

ranges   = [[-0.4,1.0],        [-0.3,0.3],        [-0.3,0.3],        [-0.5,0.5],
            [1.0,2.0],         [0.5,2.5],         [0.0,2.0],
            [-0.3,0.3],        [-0.3,0.3],        [-0.3,0.3],        [-0.5,0.5],
            [-0.5,0.5]]

vars_in  = ['u', 'v', 'w', 'T', 'p']

clipbox  = [-15.0, 10.0, -0.11, 2, -2.0, 2.0]
bbox     = [-30, 120, -2.0, 31.0, -11, 11]
# =============================================================================

dirs     = Directories( casedir )
out_dir  = dirs.pp_snp_yslice

# allocate variable that will be broadcasted to all workers

params    = None
snapdirs  = None
blocklist = None
grid3d    = None
statys    = None

# root does the preparation

if mpi.is_root:
    
    os.chdir( create_folder( out_dir ) )
    
    for var in vars_out:
        create_folder( f'./{var}/figs' )
    
    print( dirs.case_para_file)
    params = Params( dirs.case_para_file )
    
    snapdirs = sorted([dir for dir in os.listdir( dirs.snp_dir ) if os.path.isdir( os.path.join(dirs.snp_dir,dir))])
    snapdirs = [ os.path.join( dirs.snp_dir, snapdir ) for snapdir in snapdirs ]
    
    print(f"I am root, just found {len(snapdirs)} snapshots.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    
    statys = []
    stat_yslice_files = get_filelist( dirs.sup_dir, 'yslice' )
    
    for i,loc in enumerate(locs):
        staty = StatisticData( stat_yslice_files[i] )
        blocklist,_ = grid3d.select_sliced_blockgrids( 'Y', loc, bbox=bbox)
        staty.read_statistic(block_list=blocklist, vars_in=['u','v','w','p','rho'])
        statys.append(staty)
        
params    = mpi.comm.bcast( params,    root=0 )
snapdirs  = mpi.comm.bcast( snapdirs,  root=0 )
blocklist = mpi.comm.bcast( blocklist, root=0 )
grid3d    = mpi.comm.bcast( grid3d,    root=0 )
statys    = mpi.comm.bcast( statys,    root=0 )

casecode  = params.casecode
u_ref     = params.u_ref
T_ref     = params.T_ref
p_ref     = params.p_ref
delta     = params.delta_0
rho_ref   = params.rho_ref
rescale   = [-params.x_imp, 0.0, 0.0, delta, delta, delta]

mpi.barrier()

os.chdir( out_dir )
clock = timer("show slices from snapshot_Y_xxxx.bin:")

# read in the snapshots and show them

def show_slice( snapdir ):

# -- read in the snapshot and compute density gradient

    snapfiles = get_filelist( snapdir, '_Y_' )

    for i, loc in enumerate(locs):

        print(f"Processing snapshot {snapfiles[i]}")
        snap         = Snapshot( snapfiles[i] )
        snap.grid3d  = grid3d
        blocklist, _ = grid3d.select_sliced_blockgrids( 'Y', loc, bbox=bbox)
        snap.read_snapshot( block_list=blocklist, var_read=vars_in )
        snap.compute_gradients( block_list=blocklist,grads=['grad_rho'] )
        snap.compute_vars( blocklist, ['mach'] )

        # data processing block by block

        for snapbl in snap.snap_data:
            
            bl_num = snapbl.num
            statbl = statys[i].bl[statys[i].bl_nums.index(bl_num)]
            
            snapbl.df['p_fluc']   = (snapbl.df['p']  - statbl.df['p']  )/p_ref
            snapbl.df['rho_fluc'] = (snapbl.df['rho']- statbl.df['rho'])/rho_ref
            snapbl.df['u_r']      = (snapbl.df['u']  - statbl.df['u']  )/u_ref
            snapbl.df['v_r']      = (snapbl.df['v']  - statbl.df['v']  )/u_ref
            snapbl.df['w_r']      = (snapbl.df['w']  - statbl.df['w']  )/u_ref
            snapbl.df['u']        =  snapbl.df['u']  /u_ref
            snapbl.df['v']        =  snapbl.df['v']  /u_ref
            snapbl.df['w']        =  snapbl.df['w']  /u_ref
            snapbl.df['T']        =  snapbl.df['T']  /T_ref
            snapbl.df['p']        =  snapbl.df['p']  /p_ref
            snapbl.df['rho']      =  snapbl.df['rho']/rho_ref

        # prepare for the visualization

        
        dataset     = pv.MultiBlock( snap.create_vtk_multiblock(rescale=rescale) )
        dataset     = dataset.cell_data_to_point_data().combine()
        dataset     = dataset.clip_box( clipbox, invert=False )
        sepline     = dataset.contour( [0.0], scalars='u' )
        sonline     = dataset.contour( [1.0], scalars='mach' )
        sys.stdout.flush()

        # -- loop over all the output variables

        for var in vars_out:
            
            os.chdir( f'./{var}/figs' )
            
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
            p.add_mesh( sonline, color='green',  line_width=3.0 )

            p.view_vector([0.0,1.0,0.0],viewup=[0.0,0.0,1.0])
            
            p.camera.tight(view='xz')
            image = p.screenshot(return_img=True)
            p.close()

            # pass the image from pyvista to matplotlib

            image = crop_border(image)

            itstep      = snap.itstep
            itime       = snap.itime
            figname     = f"slice_Y_{itstep:08d}.png"

            fig, ax = plt.subplots(figsize=(12.8,7.2))

            img = ax.imshow(image, extent=clipbox[:2]+clipbox[4:], 
                            cmap=cmap, clim=ranges[varslist.index(var)])
            
            ax.set_ylim( [-2.0, 2.0] )
            ax.set_xlabel(r'$z/\delta_0$')
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
    
    for i, snapdir in enumerate(snapdirs):
        show_slice( snapdir )
        clock.print_progress( i, len(snapdirs) )
        
else:
    if mpi.rank == 0:
        mpi.master_distribute( snapdirs )
    else:
        while True:
            task_index = mpi.worker_receive()
            
            if task_index is None: break
            else: 
                show_slice( snapdirs[task_index] )
                clock.print_progress( task_index, len(snapdirs), rank=mpi.rank )

mpi.barrier()

if mpi.rank == 0:

    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush() 