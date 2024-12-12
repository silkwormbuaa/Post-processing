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
import pandas              as     pd
import pyvista             as     pv
import matplotlib.pyplot   as     plt
from   matplotlib.gridspec import GridSpec

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

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

# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

casedir = '/home/wencan/temp/231124'
clipbox = [-100,100, -0.11,2.0,-2,2]
colmap  = 'coolwarm'
vars    = ['u','v','w','p']

# =============================================================================

dirs = Directories( casedir )

# - enter folder: case/postprocess/dmd/

os.chdir( dirs.pp_dmd )
snap_dir = dirs.snp_dir

# - grid for packing dmd results into vtk

grid3d = GridData( dirs.grid )
grid3d.read_grid()

# - find all the blocks in the snapshot slice

snap_files    = get_filelist( snap_dir, 'snapshot_X_004.bin' )
snapshot_temp = Snapshot( snap_files[0] )
block_list    = snapshot_temp.bl_nums

print(f"Found {len(block_list)} blocks in the snapshot slice.")
sys.stdout.flush()    


# =============================================================================
# Reconstruct snapshots with selected spdmd modes:
# =============================================================================

with timer('\n - Reconstruct snapshots'):
    
    modes_temp = DMDModes()
    
    # read case parameters

    params = Params( dirs.case_para_file )
    modes_temp.case_parameters = params.params
    
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

# =============================================================================
# check the reconstructed data
# =============================================================================

# dataset = modes_temp.recons_to_vtk( vars, block_list, 'X', grid3d )

# dataset = pv.MultiBlock( dataset )
# # dataset = dataset.cell_data_to_point_data().combine()

# p = pv.Plotter(window_size=[1920,1080])

# dataset.set_active_scalars('recons_00000_u')

# p.add_mesh( dataset )

# p.show()


# =============================================================================
# check each mode
# =============================================================================

os.chdir( create_folder(dirs.pp_dmd + '/modes_figs') )

zlist     = np.linspace(-2,2,endpoint=True,num=401)
delta     = params.delta_0
casecode  = params.casecode
wallshape = define_wall_shape( zlist*delta, casecode=casecode, write=False )/delta

# - loop over each mode

for i, indx in enumerate( modes_temp.df_ind['indxes'] ):
    
    # pack 8 (by default) phases of variables into vtk dataset
    
    dataset = modes_temp.mode_to_vtk(indx, vars, block_list, 'X', grid3d, 
                                     rescale=[0,0,0,5.2,5.2,5.2] )

    dataset = pv.MultiBlock( dataset )
    dataset = dataset.cell_data_to_point_data().combine()
    dataset = dataset.clip_box( clipbox, invert=False )

    St = modes_temp.df_ind.iloc[i]['Sts']

    images     = list()
    range_vars = list()
    for var in vars:
        
        # define the data range
        varname_list = [f'mode_{indx:05d}_{var}_{j:02d}' for j in range(8)]
        ranges       = [ dataset.get_data_range(varname) for varname in varname_list ]
        range_var    = [min( [r[0] for r in ranges] ), max( [r[1] for r in ranges] )]
        range_vars.append( range_var )
        
        for j in range(8):
        
            dataset.set_active_scalars(f'mode_{indx:05d}_{var}_{j:02d}')
    
            p = pv.Plotter(off_screen=off_screen,window_size=[1920,1080],border=False)
            
            # set color map
            
            cmap = plt.get_cmap( colmap, 51 )
            cmap.set_over('red')
            cmap.set_under('blue')
            
            p.add_mesh( dataset,
                        cmap=cmap,
                        clim=range_var,
                        lighting=False,
                        show_scalar_bar=False)
            
            p.camera.tight(view='zy')
            image = p.screenshot(return_img=True)
            p.close()
            
            image = crop_border(image)
            
            images.append(image)

# -- plot mode shape in 8 phases

    images = np.array(images).reshape( len(vars), 8, *images[0].shape )
    
    for j in range(8):
        
        fig = plt.figure(figsize=(12.8,7.2))
        fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)
        gs  = GridSpec(2,2, figure=fig)
        
        for k, image in enumerate(images[:,j]):
            
            n_row = k//2
            n_col = k%2
            
            ax = fig.add_subplot(gs[n_row,n_col])
            
            img = ax.imshow(image, 
                            extent=[-2,2,-0.11,2], 
                            cmap=cmap,
                            clim=range_vars[k])

            ax.fill_between( zlist, wallshape, -0.3, color='gray' )
            ax.set_ylim( [-0.3,2])
            ax.text(1.6,1.8,vars[k])
            
            if n_row==1: ax.set_xlabel(r'$z/\delta_0$')
            if n_col==0: ax.set_ylabel(r'$y/\delta_0$')

        fig.suptitle(fr"phase {j}/8 $\pi$, St = {St.real:.4f}")

        figname = f"mode_{indx:05d}_{j:02d}.png"
        if off_screen:
            plt.savefig( figname, dpi=150)
        else:
            plt.show()
        
        plt.close()
        
        print(f"Save {figname} done.") 

# -- convert png image to gif
    
    convert_gif = f"convert -delay 100 -layers Optimize mode_{indx:05d}_*.png mode_{indx:05d}.gif"
    
    exit_code = os.system( convert_gif )
    
    if exit_code == 0: print( f"Output mode_{indx:05d} gif.\n")
    else:              print( f"Failed to output mode_{indx:05d} gif.")   
