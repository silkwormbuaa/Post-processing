#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   dmdpost_xyz.py
@Time    :   2025/06/12 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

off_screen = True

if off_screen:
    import atexit
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()
    atexit.register(vdisplay.stop)
    
import os
import sys
import time
import pickle
import numpy               as     np
import pyvista             as     pv
import matplotlib.pyplot   as     plt
from   matplotlib.gridspec import GridSpec
from   moviepy.editor      import ImageSequenceClip

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.mpi          import MPIenv
from   vista.grid         import GridData
from   vista.timer        import timer
from   vista.params       import Params
from   vista.snapshot     import Snapshot
from   vista.dmdmodes     import DMDMode, DMDModes
from   vista.directories  import Directories
from   vista.tools        import get_filelist
from   vista.tools        import crop_border
from   vista.tools        import define_wall_shape
from   vista.directories  import create_folder
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(latex=False,fontsize=15)
mpi = MPIenv()

# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

casedir  = '/home/wencan/temp/smooth_mid/'

snaptags = ['snapshot_X_004.bin',
            'snapshot_Y_000.bin',
            'snapshot_Z_001.bin']

vars_snap_smooth  = [['u','v','w','p'],['u','v','w','p'],['u','v','w','p']]
vars_snap_rough   = [['u','v','w','p'],['u','v','w','p'],['u','v','w','p']]
colmap   = 'RdBu_r'
n_frames = 16

# =============================================================================

dirs = Directories( casedir )

# - enter folder: case/postprocess/dmd/

os.chdir( dirs.pp_dmd )
snap_dir = dirs.snp_dir

# allocate variables that will be broadcasted to all works

grid3d      = None
block_lists = None
modes_temp  = None
zlist       = None
wallshape   = None
vars_snap   = None
displaces   = None
params      = None

# root does the preparation

if mpi.is_root: 

    # - grid for packing dmd results into vtk

    grid3d = GridData( dirs.grid )
    grid3d.read_grid()

    params = Params( dirs.case_para_file )
    if 'smooth' in params.casecode:
        vars_snap  = vars_snap_smooth
    else:
        vars_snap  = vars_snap_rough

    # - find all the blocks in the snapshot slice

    block_lists = []
    displaces   = [0]
    for i in range(len(snaptags)):

        snap_files    = get_filelist( snap_dir, snaptags[i] )
        snapshot_temp = Snapshot( snap_files[0] )
        snapshot_temp.read_snapshot()
        block_lists.append(snapshot_temp.bl_nums)
        
        snapshot_temp.drop_ghost()
        snapshot_temp.assemble_block()
        snapshot_len = np.array(snapshot_temp.df[vars_snap[i]].values.ravel()).size
        displaces.append( snapshot_len + displaces[-1] )
    
    # - reconstruct snapshots with selected spdmd modes:

    with timer('\n - Reconstruct snapshots'):
        
        modes_temp = DMDModes()
        
        # read modes in /dmdmodes
        
        mode_files = get_filelist( dirs.pp_dmd + '/dmdmodes' )
        
        for mode_file in mode_files:
            
            modes_temp.add_mode( DMDMode(mode_file) )
        
        # also add std dmd reconstructed data
        
        with open('reconstructed_std_dmd.pkl','rb') as f:
            
            modes_temp.recons_std_dmd = pickle.load( f )
        
        # number of modes
        
        n_modes = len(modes_temp.modes)
        
        # print selected modes properties:
        
        print(f"\nGot { n_modes } modes:")
        print([mode.indx for mode in modes_temp.modes])
        
        print(f"\nAlpha_pol of these modes:")
        print([mode.alpha_pol for mode in modes_temp.modes])
        
        print(f"\n|Alpha_pol| of these modes:")
        print([abs(mode.alpha_pol) for mode in modes_temp.modes])
        
        print(f"\nmu of these modes:")
        print([mode.mu for mode in modes_temp.modes])
        
        print(f"\n|mu| of these modes:")
        print([abs(mode.mu) for mode in modes_temp.modes])
        
        print(f"\nSt of these modes:")
        print([mode.St for mode in modes_temp.modes])
        
    # - reconstruct dmd modes and select positive modes into df_ind

        step = 1
        modes_temp.reconstruct( step )

        print(f"\nIndexes of positive modes:")
        print( modes_temp.df_ind )

    sys.stdout.flush()

    # - define the wall shape for plotting 
    
    zlist      = np.linspace(-2,2,endpoint=True,num=401)
    delta      = params.delta_0
    casecode   = params.casecode
    wallshape  = define_wall_shape( zlist*delta, casecode=casecode, write=False )/delta


grid3d     = mpi.comm.bcast( grid3d,     root=0 )
block_lists= mpi.comm.bcast( block_lists,root=0 ) 
modes_temp = mpi.comm.bcast( modes_temp, root=0 )
zlist      = mpi.comm.bcast( zlist,      root=0 )
wallshape  = mpi.comm.bcast( wallshape,  root=0 )
vars_snap  = mpi.comm.bcast( vars_snap,  root=0 )
displaces  = mpi.comm.bcast( displaces,  root=0 )
params     = mpi.comm.bcast( params,     root=0 )

clock = timer(f"dmdpost of x slices:")

os.chdir( create_folder(dirs.pp_dmd + '/modes_figs') )

# - loop over each mode

# ----------------------------------------------------------------------
# >>> postprocess_dmd_x                                          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/06/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def postprocess_dmd_x( indx ):
    
    if 'smooth' in params.casecode:
        clipbox  = [-15.0, 10.0, 0.0, 2.0, -2.0, 2.0]
    else:
        clipbox  = [-15.0, 10.0, -0.6, 4.0, -2.0, 2.0]

    df = modes_temp.df_ind.copy()
    df.reset_index( inplace=True, drop=True )
    i  = df.index[ df['indxes']==indx ].tolist()[0]
    
# -- pack 8 (by default) phases of modes to find out the data range
    
    ranges = [(float('inf'),float('-inf'))]*5
    
    for j in range(8):
        
        dataset = modes_temp.mode_to_vtk(indx, vars_snap[0], block_lists[0], 'X', grid3d,
                                         j, 8, phi_displace=displaces[0], 
                                         rescale=[0,0,0,5.2,5.2,5.2])
        dataset = pv.MultiBlock( dataset )
        
        varname_list = [f'mode_x_{indx:05d}_{var}_{j:03d}' for var in vars_snap[0]]
        range_vars   = [ dataset.get_data_range(varname) for varname in varname_list ]
        ranges       = [ ( min(ranges[k][0],range_vars[k][0]),
                           max(ranges[k][1],range_vars[k][1]) ) 
                          for k in range(len(vars_snap[0])) ]
        
# -- loop over every phase  

    St       = df.iloc[i]['Sts']
    fignames = list()

    for j in range( n_frames ):
    
        dataset = modes_temp.mode_to_vtk(indx, vars_snap[0], block_lists[0], 'X', grid3d,
                                         j, n_frames, rescale=[0,0,0,5.2,5.2,5.2])

        dataset = pv.MultiBlock( dataset )
        
        dataset = dataset.cell_data_to_point_data().combine()
        dataset = dataset.clip_box( clipbox, invert=False )

        images  = list()

        for k, var in enumerate(vars_snap[0]):
        
            dataset.set_active_scalars(f'mode_x_{indx:05d}_{var}_{j:03d}')
    
            p = pv.Plotter(off_screen=off_screen,window_size=[1920,1080],border=False)
            
            # set color map
            
            cmap = plt.get_cmap( colmap, 51 )
            cmap.set_over('red')
            cmap.set_under('blue')
            
            p.add_mesh( dataset,
                        cmap=cmap,
                        clim=ranges[k],
                        lighting=False,
                        show_scalar_bar=False)
            
            p.camera.tight(view='zy')
            image = p.screenshot(return_img=True)
            p.close()
            
            image = crop_border(image)
            
            images.append(image)

# --- plot mode shape of each variable

        fig = plt.figure(figsize=(12.8,7.2))
        fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)
        gs  = GridSpec(2,2, figure=fig)
        
        for k, image in enumerate(images):
            
            n_row = k//2
            n_col = k%2
            
            ax = fig.add_subplot(gs[n_row,n_col])
            
            img = ax.imshow(image, 
                            extent=[*clipbox[4:],*clipbox[2:4]], 
                            cmap=cmap,
                            clim=range_vars[k])

            ax.fill_between( zlist, wallshape, -0.3, color='gray' )
            ax.set_ylim( [-0.3,2])
            ax.text(1.6,1.8,vars_snap[0][k])
            
            if n_row==1: ax.set_xlabel(r'$z/\delta_0$')
            if n_col==0: ax.set_ylabel(r'$y/\delta_0$')

        fig.suptitle(fr"phase {j}/{int(n_frames/2):d} $\pi$, St = {St.real:.4f}")

        figname = f"mode_x_{indx:05d}_{j:03d}.png"
        fignames.append( figname )
        if off_screen:
            plt.savefig( figname, dpi=150)
        else:
            plt.show()
        
        plt.close()
        
        print(f"Save {figname} done.") 
        sys.stdout.flush()

# -- convert png image to mp4

    output_video = f"mode_x_{indx:05d}.mp4"
    clip = ImageSequenceClip( fignames, fps=20 )
    clip.write_videofile( output_video, codec='libx264' )
    print(f" Output video {output_video} is saved!\n")

# ----------------------------------------------------------------------
# >>> postprocess_dmd_y                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/06/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def postprocess_dmd_y( indx ):
    
    clipbox  = [-15.0, 10.0, 0.0, 4.0, -2.0, 2.0]

    df = modes_temp.df_ind.copy()
    df.reset_index( inplace=True, drop=True )
    i  = df.index[ df['indxes']==indx ].tolist()[0]
    
# -- pack 8 (by default) phases of modes to find out the data range
    
    ranges = [(float('inf'),float('-inf'))]*5
    
    for j in range(8):
        
        dataset = modes_temp.mode_to_vtk(indx, vars_snap[1], block_lists[1], 'Y', grid3d,
                                         j, 8, phi_displace=displaces[1], 
                                         rescale=[-50.4,0,0,5.2,5.2,5.2])
        dataset = pv.MultiBlock( dataset )
        
        varname_list = [f'mode_y_{indx:05d}_{var}_{j:03d}' for var in vars_snap[1]]
        range_vars   = [ dataset.get_data_range(varname) for varname in varname_list ]
        ranges       = [ ( min(ranges[k][0],range_vars[k][0]),
                           max(ranges[k][1],range_vars[k][1]) ) 
                          for k in range(len(vars_snap[1])) ]
        
# -- loop over every phase  

    St       = df.iloc[i]['Sts']
    fignames = list()

    for j in range( n_frames ):
    
        dataset = modes_temp.mode_to_vtk(indx, vars_snap[1], block_lists[1], 'Y', grid3d,
                                         j, n_frames, phi_displace=displaces[1],
                                         rescale=[-50.4,0,0,5.2,5.2,5.2])

        dataset = pv.MultiBlock( dataset )
        
        dataset = dataset.cell_data_to_point_data().combine()
        dataset = dataset.clip_box( clipbox, invert=False )

        images  = list()

        for k, var in enumerate(vars_snap[1]):
        
            dataset.set_active_scalars(f'mode_y_{indx:05d}_{var}_{j:03d}')
    
            p = pv.Plotter(off_screen=off_screen,window_size=[1920,1080],border=False)
            
            # set color map
            
            cmap = plt.get_cmap( colmap, 51 )
            cmap.set_over('red')
            cmap.set_under('blue')
            
            p.add_mesh( dataset,
                        cmap=cmap,
                        clim=ranges[k],
                        lighting=False,
                        show_scalar_bar=False)
            
            p.camera.tight(view='xz')
            image = p.screenshot(return_img=True)
            p.close()
            
            image = crop_border(image)
            
            images.append(image)

# --- plot mode shape of each variable

        fig = plt.figure(figsize=(12.8,7.2))
        fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)
        gs  = GridSpec(4,1, figure=fig)
        
        for k, image in enumerate(images):
            
            ax = fig.add_subplot(gs[k,0])
            
            img = ax.imshow(image, 
                            extent=[-15.0,10.0,-2.0,2.0], 
                            cmap=cmap,
                            clim=range_vars[k])

            ax.text(7.5,1.5,vars_snap[1][k])
            
            # do not show tick label on x axis for the first three plots
            if k<3:
                ax.set_xticklabels([])
            
            if k==3: ax.set_xlabel(r'$x/\delta_0$')
            ax.set_ylabel(r'$z/\delta_0$')

        fig.suptitle(fr"phase {j}/{int(n_frames/2):d} $\pi$, St = {St.real:.4f}")

        figname = f"mode_y_{indx:05d}_{j:03d}.png"
        fignames.append( figname )
        if off_screen:
            plt.savefig( figname, dpi=150)
        else:
            plt.show()
        
        plt.close()
        
        print(f"Save {figname} done.") 
        sys.stdout.flush()

# -- convert png image to mp4

    output_video = f"mode_y_{indx:05d}.mp4"
    clip = ImageSequenceClip( fignames, fps=20 )
    clip.write_videofile( output_video, codec='libx264' )
    print(f" Output video {output_video} is saved!\n")


# ----------------------------------------------------------------------
# >>> postprocess_dmd_z                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/06/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def postprocess_dmd_z( indx ):

    clipbox  = [-15.0, 10.0, 0.0, 6.0, -2.0, 2.0]
    
    df = modes_temp.df_ind.copy()
    df.reset_index( inplace=True, drop=True )
    i  = df.index[ df['indxes']==indx ].tolist()[0]
    
# -- pack 8 (by default) phases of modes to find out the data range
    
    ranges = [(float('inf'),float('-inf'))]*5
    
    for j in range(8):
        
        dataset = modes_temp.mode_to_vtk(indx, vars_snap[2], block_lists[2], 'Z', grid3d,
                                         j, 8, phi_displace=displaces[2], 
                                         rescale=[-50.4,0,0,5.2,5.2,5.2])
        dataset = pv.MultiBlock( dataset )
        
        varname_list = [f'mode_z_{indx:05d}_{var}_{j:03d}' for var in vars_snap[2]]
        range_vars   = [ dataset.get_data_range(varname) for varname in varname_list ]
        ranges       = [ ( min(ranges[k][0],range_vars[k][0]),
                           max(ranges[k][1],range_vars[k][1]) ) 
                          for k in range(len(vars_snap[2])) ]
        
# -- loop over every phase  

    St       = df.iloc[i]['Sts']
    fignames = list()

    for j in range( n_frames ):
    
        dataset = modes_temp.mode_to_vtk(indx, vars_snap[2], block_lists[2], 'Z', grid3d,
                                         j, n_frames, phi_displace=displaces[2],
                                         rescale=[-50.4,0,0,5.2,5.2,5.2])
        dataset = pv.MultiBlock( dataset )
        
        dataset = dataset.cell_data_to_point_data().combine()
        dataset = dataset.clip_box( clipbox, invert=False )

        images  = list()

        for k, var in enumerate(vars_snap[2]):
        
            dataset.set_active_scalars(f'mode_z_{indx:05d}_{var}_{j:03d}')
    
            p = pv.Plotter(off_screen=off_screen,window_size=[1920,1080],border=False)
            
            # set color map
            
            cmap = plt.get_cmap( colmap, 51 )
            cmap.set_over('red')
            cmap.set_under('blue')
            
            p.add_mesh( dataset,
                        cmap=cmap,
                        clim=ranges[k],
                        lighting=False,
                        show_scalar_bar=False)
            
            p.camera.tight(view='xy')
            image = p.screenshot(return_img=True)
            image = crop_border(image)
            images.append(image)

            p.close()
            
# --- plot mode shape of each variable

        fig = plt.figure(figsize=(12.8,7.2))
        fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)
        gs  = GridSpec(2,2, figure=fig)
        
        for k, image in enumerate(images):
            
            n_row = k//2
            n_col = k%2
            
            ax = fig.add_subplot(gs[n_row,n_col])
            
            img = ax.imshow(image, 
                            extent=[-15,10,0,6], 
                            cmap=cmap,
                            clim=range_vars[k])

            ax.text(7.5,4,vars_snap[2][k])
            
            if n_row==1: ax.set_xlabel(r'$x/\delta_0$')
            if n_col==0: ax.set_ylabel(r'$y/\delta_0$')

        fig.suptitle(fr"phase {j}/{int(n_frames/2):d} $\pi$, St = {St.real:.4f}")

        figname = f"mode_z_{indx:05d}_{j:03d}.png"
        fignames.append( figname )
        if off_screen:
            plt.savefig( figname, dpi=150)
        else:
            plt.show()
        
        plt.close()
        
        print(f"Save {figname} done.") 
        sys.stdout.flush()

# -- convert png image to mp4

    output_video = f"mode_z_{indx:05d}.mp4"
    clip = ImageSequenceClip( fignames, fps=20 )
    clip.write_videofile( output_video, codec='libx264' )
    print(f" Output video {output_video} is saved!\n")
    
    
# =============================================================================
# do post processing in parallel
# =============================================================================

if mpi.size == 1:
    
    print("No worker available. Master should do all tasks.")
    
    for i, indx in enumerate( modes_temp.df_ind['indxes'][1:] ):
        
        postprocess_dmd_x( indx )
        postprocess_dmd_y( indx )
        postprocess_dmd_z( indx )
        clock.print_progress( i, len(modes_temp.df_ind['indxes']))
        
else:
    if mpi.is_root:
        
        indxes = np.array( modes_temp.df_ind['indxes'] )
        mpi.master_distribute( indxes )
    
    else:
        while True:
            task_index = mpi.worker_receive()
            
            if task_index is None: break
            else:
                indxes = np.array( modes_temp.df_ind['indxes'] )
                postprocess_dmd_x( indxes[task_index] )
                postprocess_dmd_y( indxes[task_index] )
                postprocess_dmd_z( indxes[task_index] )
                clock.print_progress( task_index, 
                                      len(modes_temp.df_ind['indxes']),
                                      rank=mpi.rank)

mpi.barrier()

if mpi.is_root:
    
    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush() 
