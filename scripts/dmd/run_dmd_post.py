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

import numpy             as     np

import pandas            as     pd

import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.timer       import timer

from   vista.snapshot    import Snapshot

from   vista.dmdmodes    import DMDMode, DMDModes

from   scipy.interpolate import griddata

from   vista.tools       import get_filelist

from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_dmd_mode 


# read in one snapshot file to get grid vectors

snap_dir = os.getcwd()

step = 1

with timer('\nGet snapshots file and grid vector'):
    
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
    
    

# match data and grids:

with timer('\nReconstruct snapshots'):
    
    mode_files = get_filelist( snap_dir + '/dmdmodes' )
    
    modes_temp = DMDModes()
    
    for mode_file in mode_files:
        
        modes_temp.add_mode( DMDMode(mode_file) )
        
    # print selected modes properties:
    
    print(f"\nGot these modes:")
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
    


# read case parameters

modes_temp.case_parameters = read_case_parameter('case_parameters')


with timer("\nShow snapshots"):
    
    modes_temp.match_mesh( df, snap_type )
    
    print(modes_temp.df_modes)
    
    xmin = -50.0 #modes_temp.df_modes['x'].min()
    xmax = 100.0 #modes_temp.df_modes['x'].max()
    ymin = modes_temp.df_modes['y'].min()
    ymax = 40.0  #modes_temp.df_modes['y'].max()
    
    x_i = np.linspace( xmin, xmax, 451)
    y_i = np.linspace( ymin, ymax, 121)
    
    modes_temp.grids_interp = np.meshgrid( x_i, y_i )
    
    
    # plot reconstructed data
    
    xx,yy,v = modes_temp.interp_recons( 0 )

    plot_dmd_mode( (xx,yy), v )
    
    
    # plot dmd modes
    
    n_mode = 0
    
    xx,yy,v = modes_temp.interp_mode( n_mode )
    
    plot_dmd_mode( (xx,yy), v )
    

    
    
    
    
        

