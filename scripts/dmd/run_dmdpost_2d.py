#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''
@File    :   run_dmd_post.py
@Time    :   2023/06/05 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import time
import cmath
import pickle
import numpy             as     np
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.dmdmodes    import DMDMode, DMDModes
from   vista.tools       import get_filelist
from   vista.plot_style  import plot_dmd_mode 
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

snap_dir = os.getcwd()

step = 1

t_0 = time.time()

with timer('\n - Get snapshots file and grid vector'):
    
    snap_files = get_filelist( snap_dir + '/snapshots' )
    
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

    df = pd.DataFrame( GX, columns = GX_header )

sys.stdout.flush()    
    
# =============================================================================
# Reconstruct snapshots with selected spdmd modes:
# =============================================================================

with timer('\n - Reconstruct snapshots'):
    
    modes_temp = DMDModes()
    
    # read case parameters

    params = Params( 'case_parameters' )
    modes_temp.case_parameters = params.params
    
    # read modes in /dmdmodes
    
    mode_files = get_filelist( snap_dir + '/dmdmodes' )
    
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
# match mesh and reconstructed data( both std_dmd and spdmd ), each modes
# =============================================================================
    
# match mesh with data
with timer("\n - Match mesh "):

    modes_temp.match_mesh( df, snap_type )

    print(modes_temp.df_modes)

sys.stdout.flush()

# =============================================================================
# generate the grid that will be interpolated on
# =============================================================================

with timer("\n - Generate the interpolation grid "):

    if snap_type == 'Z':

        x1min = -50.0 #modes_temp.df_modes['x'].min()
        x1max = 100.0 #modes_temp.df_modes['x'].max()
        x2min = modes_temp.df_modes['y'].min()
        x2max = 40.0  #modes_temp.df_modes['y'].max()
        
    elif snap_type == 'W' or snap_type == 'Y':
        
        x1min = modes_temp.df_modes[GX_header[0]].min()
        x1max = modes_temp.df_modes[GX_header[0]].max()
        x2min = modes_temp.df_modes[GX_header[1]].min()
        x2max = modes_temp.df_modes[GX_header[1]].max()

    else:
        
        raise ValueError("snap_type not supported yet.")
        
    x_1 = np.linspace( x1min, x1max, 451)
    x_2 = np.linspace( x2min, x2max, 121)

    modes_temp.grids_interp = np.meshgrid( x_1, x_2 )

sys.stdout.flush()

# =============================================================================
# interpolate and output figures
# =============================================================================

with timer("\n - Interpolate and output figures "):

    # enter ./modesfigures directory

    if os.path.exists("./modesfigures"): os.chdir("./modesfigures")
        
    else: os.mkdir("./modesfigures"); os.chdir("./modesfigures")


    # plot reconstructed flow field with selected spdmd modes

    header = "recons_00000"

    xx,yy,v = modes_temp.interp_recons( header, snap_type )

    plot_dmd_mode( (xx,yy), v, 
                    filename="reconstructed_spdmd.png", 
                    colorbar=True )


    # plot reconstructed flow field with std dmd modes

    header = "recons_std_dmd"

    xx,yy,v = modes_temp.interp_recons( header, snap_type )

    plot_dmd_mode( (xx,yy), v, 
                    filename="reconstructed_std_dmd.png", 
                    colorbar=True )


    # indexes of positive modes

    indxes = np.array( modes_temp.df_ind['indxes'] )

    print(indxes)


    # plot mean mode

    print(f"The index of max alpha_pol is {indxes[0]}.\n")

    xx,yy,v = modes_temp.interp_mode( indxes[0] )

    plot_dmd_mode( (xx,yy), v,
                    filename="mode_mean.png",
                    colorbar=True,
                    title="mean mode")


    # plot oscillating dmd modes
        
    phases = [cmath.rect(1.0, cmath.pi*0.25*i) for i in range(8)]
    
    Sts = np.array( modes_temp.df_ind['Sts'] )

    for n in range( len(indxes) ):
        
        for i, phase in enumerate(phases):
            
            xx,yy,v = modes_temp.interp_mode( indxes[n], snap_type, phase=phase)
            
            filename = f"mode_{n:02d}_{indxes[n]:05d}_{i}"
            
            title = f"St={Sts[n]:.3f} "+ f"phase = {i}/4"+r"$\pi$"
            
            if i == 0:
                
                cmax = np.max( np.abs(v) )*1.3
                
                clevel = np.linspace( -cmax, cmax, 51)
            
            plot_dmd_mode( (xx,yy), v,
                            filename=filename+'.png',
                            colorbar=True,
                            clevel=clevel,
                            title=title)

        # convert png image to gif
        
        convert_gif = f"convert -delay 100 mode_{n:02d}*.png mode_{n:02d}.gif"
        
        exit_code = os.system( convert_gif )
        
        if exit_code == 0: print( f"Output mode_{n:02d} gif.")
        else:              print( f"Failed to output mode_{n:02d} gif.")
        
                
t_end = time.time()

print(f"\n - dmdpost totally took {t_end-t_0:8.2f}s.")

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()       
