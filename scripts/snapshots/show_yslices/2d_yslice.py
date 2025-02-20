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

casedir  = '/home/wencan/temp/231124/'

locs     = [0.0]
locs     = np.array(locs) * 5.2 + 50.4

vars_out = ['u',                'v',               'w',                'p',
            'T',                'rho',             'mach',             'w1',
            'u_r',              'v_r',             'w_r',              'p_fluc',
            'rho_fluc',         'DS',              'w1_r']

varslist = ['u',                'v',               'w',                'p',
            'T',                'rho',             'mach',             'w1',
            'u_r',              'v_r',             'w_r',              'p_fluc',
            'rho_fluc',         'DS',              'w1_r']

labels   = [r'$u/u_{\infty}$',        r'$v/u_{\infty}$',       r'$w/u_{\infty}$',     r'$p/p_{\infty}$',
            r'$T/T_{\infty}$',        r'$\rho/\rho_{\infty}$', r'$Mach$',             r'$\omega$',
            r'$u_{r}/u_{\infty}$',    r'$v_{r}/u_{\infty}$',   r'$w_{r}/u_{\infty}$', r"$p'/p_{\infty}$",
            r"$\rho'/\rho_{\infty}$", r"$DS$",                 r"$\omega_r$"]

colmaps  = ['coolwarm',        'coolwarm',        'coolwarm',        'coolwarm',
            'plasma',          'plasma',          'coolwarm',        'coolwarm',
            'coolwarm',        'coolwarm',        'coolwarm',        'coolwarm',
            'coolwarm',        'Greys_r',         'coolwarm']

ranges   = [[-0.4,1.0],        [-0.3,0.3],        [-0.3,0.3],        [-0.5,0.5],
            [1.0,2.0],         [0.5,2.5],         [0.0,2.0],         [-6,6],
            [-0.3,0.3],        [-0.3,0.3],        [-0.3,0.3],        [-0.5,0.5],
            [-0.5,0.5],        [0.0,0.8],         [-6,6]]

vars_in  = ['u', 'v', 'w', 'T', 'p']

clipbox  = [-999, 999, -0.11, 2, -2.0, 2.0]
bbox     = [-50, 120, -2.0, 31.0, -11, 11]
# =============================================================================

dirs     = Directories( casedir )
out_dir  = dirs.pp_snp_xslice

# allocate variable that will be broadcasted to all workers

params    = None
snapdirs  = None
grid3d    = None
statxs    = None

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
    
    statxs = []
    stat_xslice_files = get_filelist( dirs.sup_dir, 'xslice' )
    
    for i,loc in enumerate(locs):
        statx = StatisticData( stat_xslice_files[i] )
        blocklist,_ = grid3d.select_sliced_blockgrids( 'X', loc, bbox=bbox)
        statx.read_statistic(block_list=blocklist, vars_in=['u','v','w','p','rho'])
        statxs.append(statx)
        
params    = mpi.comm.bcast( params,    root=0 )
snapdirs  = mpi.comm.bcast( snapdirs,  root=0 )
blocklist = mpi.comm.bcast( blocklist, root=0 )
grid3d    = mpi.comm.bcast( grid3d,    root=0 )
statxs    = mpi.comm.bcast( statxs,    root=0 )

casecode  = params.casecode
u_ref     = params.u_ref
T_ref     = params.T_ref
p_ref     = params.p_ref
delta     = params.delta_0
rho_ref   = params.rho_ref
rescale   = [-params.x_imp, 0.0, 0.0, delta, delta, delta]

mpi.barrier()

os.chdir( out_dir )
clock = timer("show slices from snapshot_Z_xxxx.bin:")

zlist = np.linspace(-2,2,endpoint=True,num=401)
wallshape = define_wall_shape( zlist*delta, casecode=casecode, write=False )/delta

# read in the snapshots and show them

def show_slice( snapdir ):

# -- read in the snapshot and compute density gradient

    images = []
    datasets = []
    seplines = []
    sonlines = []

    snapfiles = get_filelist( snapdir, '_X_' )

    for i, loc in enumerate(locs):

        print(f"Processing snapshot {snapfiles[i]}")
        snap         = Snapshot( snapfiles[i] )
        snap.grid3d  = grid3d
        blocklist, _ = grid3d.select_sliced_blockgrids( 'X', loc, bbox=bbox)
        snap.read_snapshot( block_list=blocklist, var_read=vars_in )
        snap.compute_gradients( block_list=blocklist,grads=['grad_rho','vorticity'] )
        snap.compute_vars( blocklist, ['mach'] )

        # data processing block by block

        for snapbl in snap.snap_data:
            
            bl_num = snapbl.num
            statbl = statxs[i].bl[statxs[i].bl_nums.index(bl_num)]
            
            snapbl.df['p_fluc']   = (snapbl.df['p']  - statbl.df['p']  )/p_ref
            snapbl.df['rho_fluc'] = (snapbl.df['rho']- statbl.df['rho'])/rho_ref
            snapbl.df['u_r']      = (snapbl.df['u']  - statbl.df['u']  )/u_ref
            snapbl.df['v_r']      = (snapbl.df['v']  - statbl.df['v']  )/u_ref
            snapbl.df['w_r']      = (snapbl.df['w']  - statbl.df['w']  )/u_ref
            snapbl.df['u']        =  snapbl.df['u']  /u_ref
            snapbl.df['v']        =  snapbl.df['v']  /u_ref
            snapbl.df['w']        =  snapbl.df['w']  /u_ref
            snapbl.df['w1']       =  snapbl.df['w1'] /u_ref*delta
            snapbl.df['T']        =  snapbl.df['T']  /T_ref
            snapbl.df['p']        =  snapbl.df['p']  /p_ref
            snapbl.df['rho']      =  snapbl.df['rho']/rho_ref
            snapbl.df['DS']       =  compute_DS( snapbl.df['grad_rho'], min=0.0, max=2.0)
            
            # compute the relative vorticity

            g   = grid3d.g[bl_num-1] 

            v_r = np.array(snapbl.df['v_r']).reshape( snapbl.npz, snapbl.npy )
            w_r = np.array(snapbl.df['w_r']).reshape( snapbl.npz, snapbl.npy )
            
            dwr_dy = np.gradient( w_r, g.gy, axis=0 )
            dvr_dz = np.gradient( v_r, g.gz, axis=1 )
            
            w1_r = (dwr_dy - dvr_dz)*delta
            
            snapbl.df['w1_r'] = w1_r.flatten()
            
        # prepare for the visualization

        
        dataset     = pv.MultiBlock( snap.create_vtk_multiblock(rescale=rescale) )
        dataset     = dataset.cell_data_to_point_data().combine()
        dataset     = dataset.clip_box( clipbox, invert=False )
        sepline     = dataset.contour( [0.0], scalars='u' )
        sonline     = dataset.contour( [1.0], scalars='mach' )
        datasets.append( dataset )
        seplines.append( sepline )
        sonlines.append( sonline )
        sys.stdout.flush()

    # -- loop over all the output variables

        for var in vars_out:
            
            os.chdir( f'./{var}/figs' )
            
            for i, dataset in enumerate(datasets):
            
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

                p.view_vector([1.0,0.0,0.0],viewup=[0.0,1.0,0.0])
                
                p.camera.tight(view='zy')
                image = p.screenshot(return_img=True)
                p.close()

                # pass the image from pyvista to matplotlib

                image = crop_border(image)
                images.append( image )

            itstep      = snap.itstep
            itime       = snap.itime
            figname     = f"slice_X_{itstep:08d}.png"

            fig, ax = plt.subplots(figsize=(12.8,7.2))

            img = ax.imshow(image, extent=clipbox[4:]+clipbox[2:4], 
                            cmap=cmap, clim=ranges[varslist.index(var)])
            
            ax.fill_between( zlist, wallshape, -0.3, color='gray')
            ax.set_ylim( [-0.3, clipbox[3]] )
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