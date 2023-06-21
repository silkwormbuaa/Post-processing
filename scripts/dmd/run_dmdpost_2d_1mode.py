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

import cmath

import time

import numpy             as     np

import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.timer       import timer

from   vista.snapshot    import Snapshot

from   vista.dmdmodes    import DMDMode, DMDModes

from   vista.tools       import get_filelist

from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_dmd_mode 

t_0 = time.time()
# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

snap_dir = os.getcwd()

step = 1

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
    
    
# =============================================================================
# Reconstruct snapshots with selected spdmd modes:
# =============================================================================

with timer('\n - Reconstruct snapshots'):
    
    modes_temp = DMDModes()
    
    # read case parameters

    modes_temp.case_parameters = read_case_parameter('case_parameters')
    
    # read modes in /dmdmodes
    
    mode_files = get_filelist( snap_dir + '/dmdmodes' )
    
    for mode_file in mode_files:
        
        modes_temp.add_mode( DMDMode(mode_file) )
    
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

# =============================================================================
# match mesh and reconstructed data( both std_dmd and spdmd ), each modes
# =============================================================================
    
# match mesh with data
with timer("\n - Match mesh "):

    modes_temp.match_mesh( df, snap_type )

    print(modes_temp.df_modes)



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


# =============================================================================
# interpolate and output figures
# =============================================================================

with timer("\n - Interpolate and output figures "):

    # enter ./modesfigures directory

    if os.path.exists("./modesfigures"): os.chdir("./modesfigures")
        
    else: os.mkdir("./modesfigures"); os.chdir("./modesfigures")


    # indexes of positive modes

    indxes = np.array( modes_temp.df_ind['indxes'] )

    print(indxes)


    # plot oscillating dmd modes
        
    phases = [cmath.rect(1.0, cmath.pi*0.25*i) for i in range(8)]

    n = 3
        
    for i, phase in enumerate(phases):
        
        xx,yy,v = modes_temp.interp_mode( indxes[n], phase=phase )
        
        filename = f"mode_{n:02d}_{indxes[n]:05d}-unicbar_{i}"
        
        title = f"phase = {i}/4"+r"$\pi$"
        
        if i == 0:
            
            cmax = np.max( np.abs(v) )*1.1
            
            clevel = np.linspace( -0.16, 0.16, 51)
        
        plot_dmd_mode( (xx,yy), v,
                        filename=filename+'.png',
                        colorbar=True,
                        clevel=clevel,
                        title=title)

t_end = time.time()

print(f"\n - dmdpost totally took {t_end-t_0:8.2f}s.")