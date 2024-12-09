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
import cmath
import pickle
import pandas             as     pd
import pyvista            as     pv
import matplotlib.pyplot  as     plt

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
from   vista.directories  import create_folder
from   vista.plot_setting import set_plt_rcparams
#from   vista.log         import Logger
#sys.stdout = Logger( os.path.basename(__file__) )
set_plt_rcparams(latex=False,fontsize=15)

# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

casedir = '/home/wencan/temp/231124'
clipbox = [-100,100, -0.3,2.0,-2,2]
colmap  = 'coolwarm'

# =============================================================================


dirs = Directories( casedir )

os.chdir( dirs.pp_dmd )

snap_dir = dirs.snp_dir

grid3d = GridData( dirs.grid )
grid3d.read_grid()

step = 1

t_0 = time.time()

with timer('\n - Get snapshots file and grid vector'):
    
    snap_files = get_filelist( snap_dir, 'snapshot_X_004.bin' )
    
    testfile = snap_files[0]
    
    snapshot_temp = Snapshot( testfile )
    
    if snapshot_temp.type == 'slice': 
        
        snap_type = snapshot_temp.slic_type
        if   snap_type == 'X': GX_header=['y','z']
        elif snap_type == 'Z': GX_header=['x','y']
        elif snap_type == 'W' or snap_type == 'Y': GX_header=['x','z']      
        
    elif snapshot_temp.type == 'block': 
        
        snap_type = 'block'
        GX_header = ['x','y','z']
        
    
    GX = snapshot_temp.get_grid_vectors( buff=3 )
    
    block_list = snapshot_temp.bl_nums
    
    df = pd.DataFrame( GX, columns = GX_header )

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
    
    modes_temp.reconstruct( step )

    print(f"\nIndexes of positive modes:")
    print( modes_temp.df_ind )

sys.stdout.flush()

# =============================================================================
# check the reconstructed data
# =============================================================================

# dataset = modes_temp.recons_to_vtk( ['u','v','w','p'], block_list, 'X', grid3d )

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


for i, indx in enumerate( modes_temp.df_ind['indxes'] ):
    
    dataset = modes_temp.mode_to_vtk(indx, ['u','v','w','p'], block_list, 'X', grid3d, 
                                     rescale=[0,0,0,5.2,5.2,5.2] )

    dataset = pv.MultiBlock( dataset )
    dataset = dataset.cell_data_to_point_data().combine()
    dataset = dataset.clip_box( clipbox, invert=False )

    St = modes_temp.df_ind.iloc[i]['Sts']

    for var in ['u','v','w','p']:
        
        # define the data range
        varname_list = [f'mode_{indx:05d}_{j:02d}_{var}' for j in range(8)]
        ranges       = [ dataset.get_data_range(varname) for varname in varname_list ]
        range_var    = [min( [r[0] for r in ranges] ), max( [r[1] for r in ranges] )]
        
        for j in range(8):
        
            dataset.set_active_scalars(f'mode_{indx:05d}_{j:02d}_{var}')
    
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
            
            fig,ax = plt.subplots(figsize=(12.8,7.2))
            
            img = ax.imshow(image, 
                            extent=[-2,2,-0.3,2], 
                            cmap='coolwarm',
                            clim=range_var)

            ax.set_xlabel(r'$z/\delta_0$')
            ax.set_ylabel(r'$y/\delta_0$')

            cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5 ) 

            plt.title(fr"phase {j}/8 $\pi$, St = {St:8.4f}", loc='center',y=1.05)


            figname = f"mode_{indx:05d}_{j:02d}_{var}.png"
            if off_screen:
                plt.savefig( figname, dpi=150)
            else:
                plt.show()
            
            plt.close()
            
            print(f"Save {figname} done.")      
            


