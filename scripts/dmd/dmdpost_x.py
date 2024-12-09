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


import os
import sys
import time
import cmath
import pickle
import numpy             as     np
import pandas            as     pd
import pyvista           as     pv

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.dmdmodes    import DMDMode, DMDModes
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.plot_style  import plot_dmd_mode 
from   vista.directories import create_folder
#from   vista.log         import Logger
#sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

casedir = '/home/wencan/temp/231124'

dirs = Directories( casedir )

os.chdir( dirs.pp_dmd)

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
# match mesh and reconstructed data( both std_dmd and spdmd ), each modes
# =============================================================================

dataset = modes_temp.recons_to_vtk( ['u','v','w','p'], block_list, 'X', grid3d )



dataset = pv.MultiBlock( dataset )
# dataset = dataset.cell_data_to_point_data().combine()

p = pv.Plotter(window_size=[1920,1080])

dataset.set_active_scalars('recons_00000_u')

p.add_mesh( dataset )

p.show()

