#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   dmdpost_x.py
@Time    :   2024/12/09 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

off_screen = True

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()

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
from   vista.tools        import crop_to_rect_map
from   vista.directories  import create_folder
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(latex=False,fontsize=15)
mpi = MPIenv()

# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

casedir  = '/home/wencan/temp/231124/'
clipbox  = [-100.0, 150.0, 0.0, 2.0, -2.0, 2.0]
colmap   = 'coolwarm'
vars     = ['u','v','w','p']
n_frames = 64
snaptag  = 'snapshot_X_004.bin'

# =============================================================================

dirs = Directories( casedir )

# - enter folder: case/postprocess/dmd/

os.chdir( dirs.pp_dmd )
snap_dir = dirs.snp_dir

# allocate variables that will be broadcasted to all works

grid3d     = None
block_list = None
modes_temp = None
zlist      = None
wallshape  = None

# root does the preparation

if mpi.is_root: 

    # - grid for packing dmd results into vtk

    grid3d = GridData( dirs.grid )
    grid3d.read_grid()

    # - find all the blocks in the snapshot slice

    snap_files    = get_filelist( snap_dir, snaptag )
    snapshot_temp = Snapshot( snap_files[0] )
    snapshot_temp.read_snapshot()
    block_list    = snapshot_temp.bl_nums

    print(f"Found {len(block_list)} blocks in the snapshot slice.")
    sys.stdout.flush()    

    # - reconstruct snapshots with selected spdmd modes:

    with timer('\n - Reconstruct snapshots'):
        
        modes_temp = DMDModes()
        
        # read case parameters

        params = Params( dirs.case_para_file )
        
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
    clipbox[2] = params.y_min 
    delta      = params.delta_0
    casecode   = params.casecode
    wallshape  = define_wall_shape( zlist*delta, casecode=casecode, write=False )/delta


grid3d     = mpi.comm.bcast( grid3d,     root=0 )
block_list = mpi.comm.bcast( block_list, root=0 ) 
modes_temp = mpi.comm.bcast( modes_temp, root=0 )
clipbox    = mpi.comm.bcast( clipbox,    root=0 )
zlist      = mpi.comm.bcast( zlist,      root=0 )
wallshape  = mpi.comm.bcast( wallshape,  root=0 )


clock = timer(f"dmdpost of x slices:")

os.chdir( create_folder(dirs.pp_dmd + '/modes_figs') )

# - loop over each mode

def postprocess_dmd_x( indx ):

    df = modes_temp.df_ind.copy()
    df.reset_index( inplace=True, drop=True )
    i  = df.index[ df['indxes']==indx ].tolist()[0]
    
# -- pack 8 (by default) phases of modes to find out the data range
    
    ranges = [(float('inf'),float('-inf'))]*5
    
    for j in range(8):
        
        dataset = modes_temp.mode_to_vtk(indx, vars, block_list, 'X', grid3d,
                                         j, 8, rescale=[0,0,0,5.2,5.2,5.2])
        dataset = pv.MultiBlock( dataset )
        
        varname_list = [f'mode_{indx:05d}_{var}_{j:03d}' for var in vars]
        range_vars   = [ dataset.get_data_range(varname) for varname in varname_list ]
        ranges       = [ ( min(ranges[k][0],range_vars[k][0]),
                           max(ranges[k][1],range_vars[k][1]) ) 
                          for k in range(len(vars)) ]
        
# -- loop over every phase  

    St       = df.iloc[i]['Sts']
    fignames = list()

    for j in range( n_frames ):
    
        dataset = modes_temp.mode_to_vtk(indx, vars, block_list, 'X', grid3d,
                                         j, n_frames, rescale=[0,0,0,5.2,5.2,5.2])

        dataset = pv.MultiBlock( dataset )
        
        
        dataset = dataset.cell_data_to_point_data().combine()
        dataset = dataset.clip_box( clipbox, invert=False )
        
        # set vector [u,v,w]
        
        vel_vars = [f'mode_{indx:05d}_{var}_{j:03d}' for var in ['u','v','w']]
        
        dataset = add_velocity_vector( vel_vars, dataset )
        
#        sys.exit()
        
        # streamlines = dataset.streamlines(
        #      pointa=(0,1.0,-2),
        #      pointb=(0,1.0, 2),
        #      n_points=25,
        #      max_time=500,
        #      integration_direction='both',
        #      surface_streamlines=True,
        # )

        images  = list()

        for k, var in enumerate(vars):
        
            dataset.set_active_scalars(f'mode_{indx:05d}_{var}_{j:03d}')
    
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
            
            if k == 0:
                
                p.add_arrows( dataset.points, dataset['velocity'],
                              mag=40, color='black')
                # p.add_mesh(streamlines.tube(radius=0.01), 
                #         color='black',
                #         line_width=0.01)
            
            p.camera.tight(view='zy')
            image = p.screenshot(return_img=True)
            p.close()
            
            image = crop_to_rect_map(image, buff=500)
            
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
            ax.text(1.6,1.8,vars[k])
            
            if n_row==1: ax.set_xlabel(r'$z/\delta_0$')
            if n_col==0: ax.set_ylabel(r'$y/\delta_0$')

        fig.suptitle(fr"phase {j}/{int(n_frames/2):d} $\pi$, St = {St.real:.4f}")

        figname = f"mode_{indx:05d}_{j:03d}.png"
        fignames.append( figname )
        if off_screen:
            plt.savefig( figname, dpi=150)
        else:
            plt.show()
        
        plt.close()
        
        print(f"Save {figname} done.") 
        sys.stdout.flush()

# -- convert png image to mp4

    output_video = f"mode_{indx:05d}.mp4"
    clip = ImageSequenceClip( fignames, fps=20 )
    clip.write_videofile( output_video, codec='libx264' )
    print(f" Output video {output_video} is saved!\n")    


# -- convert png image to gif
    
    # convert_gif = f"convert -delay 10 -layers Optimize mode_{indx:05d}_*.png mode_{indx:05d}.gif"
    
    # exit_code = os.system( convert_gif )
    
    # if exit_code == 0: print( f"Output mode_{indx:05d} gif.\n")
    # else:              print( f"Failed to output mode_{indx:05d} gif.")   


def add_velocity_vector( vars, dataset:pv.MultiBlock ) -> pv.MultiBlock:
    
    """
    add [u,v,w] as vector 'velocity' for pv.MultiBlock
    """
    
    u = np.zeros_like(dataset[vars[0]])
    v = dataset[vars[1]]
    w = dataset[vars[2]]
    
    vector = np.array([u,v,w]).T

    mask         = np.zeros_like( vector, dtype=bool )
    mask[::10,::] = True
    vector       = vector*mask
        
    dataset['velocity'] = vector
    
    return dataset


# =============================================================================
# do post processing in parallel
# =============================================================================

if mpi.size == 1:
    
    print("No worker available. Master should do all tasks.")
    
    for i, indx in enumerate( modes_temp.df_ind['indxes'][1:] ):
        
        postprocess_dmd_x( indx )
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
                clock.print_progress( task_index, 
                                      len(modes_temp.df_ind['indxes']),
                                      rank=mpi.rank)

mpi.barrier()

if mpi.is_root:
    
    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush() 
