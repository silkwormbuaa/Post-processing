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
        
    
    print(f"\nGot these modes:")
    print([mode.indx for mode in modes_temp.modes])
    
    
    modes_temp.reconstruct( step )
    
    
with timer("\nShow snapshots"):
    
    modes_temp.match_mesh( df, snap_type )
    
    print(modes_temp.df_modes)
    
    xmin = modes_temp.df_modes['x'].min()
    xmax = modes_temp.df_modes['x'].max()
    ymin = modes_temp.df_modes['y'].min()
    ymax = modes_temp.df_modes['y'].max()
    
    x_i = np.linspace( xmin, xmax, 513)
    y_i = np.linspace( ymin, ymax, 129)
    
    xx, yy = np.meshgrid( x_i, y_i )
    
    dataname = [f"recons_{i:05d}" for i in range(step)]

    
    # plot reconstructed data
    
    for i in range( step ):
        
        p = griddata((modes_temp.df_modes['x'],modes_temp.df_modes['y']),
                    np.array(modes_temp.df_modes[dataname[i]]).real,
                    (xx,yy), method='linear')
        
        fig, ax = plt.subplots()
        
    #    p = np.array(modes_temp.df_modes['recons_00000']).reshape((128,512))
        
        plt.imshow(p.real,extent=(xmin,xmax,ymin,ymax),origin='lower')
        
        ax.set_title('reconstructed mean')
        
        plt.show()
    
    
    # plot dmd modes
    
    dataname = [f"phi_{i:05d}" for i in modes_temp.indxes]
    
    for i in range(5):
        
        phi = griddata((modes_temp.df_modes['x'],modes_temp.df_modes['y']),
                    np.array(modes_temp.df_modes[dataname[i]]).real,
                    (xx,yy), method='linear')
        
        fig, ax = plt.subplots()
        
    #    p = np.array(modes_temp.df_modes['recons_00000']).reshape((128,512))
        
        plt.imshow(phi.real,extent=(xmin,xmax,ymin,ymax),origin='lower')
        
        ax.set_title(f'phi{i:05d}')
        
        plt.show()
    
    
    
    
        

