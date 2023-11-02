#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   run_dmdpost_2d_combined.py
@Time    :   2023/08/25 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import cmath
import time
import pickle
import copy

import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.timer       import timer

from   vista.snapshot    import Snapshot

from   vista.dmdmodes    import DMDMode, DMDModes

from   vista.tools       import get_filelist
from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_dmd_mode 
from   vista.plot_style  import plot_combined_dmd_mode

from   vista.log         import Logger
sys.stdout = Logger(filename=os.path.basename(__file__))


# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

snap_dir = os.getcwd()

step = 1

t_0 = time.time()

with timer('\n - Get snapshots file and grid vector'):

# -- grid of Y

    snap_files = get_filelist( snap_dir + '/snapshots','snapshot_Y_003.bin' )
    
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
    dfY = pd.DataFrame( GX, columns = GX_header )

# -- grid of Z 

    snap_files = get_filelist( snap_dir + '/snapshots','snapshot_Z_001.bin' )
    
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
    dfZ = pd.DataFrame( GX, columns = GX_header )

sys.stdout.flush() 

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

    # also add std dmd reconstructed data
    
    with open('reconstructed_std_dmd.pkl','rb') as f:
        
        modes_temp.recons_std_dmd = pickle.load( f )
    
    # number of modes
    
    n_modes = len(modes_temp.modes)

# -- read in the length of snapshots on each rank and distribute Phis to 
#    where they belong to (Y, or Z 's Phis). 
    
    len_snapshots_file = snap_dir + '/len_snapshots.pkl'
    if os.path.exists( len_snapshots_file ):
        with open(len_snapshots_file, 'rb' ) as f:
            len_Y = pickle.load( f )
            len_Z = pickle.load( f )
    else: raise FileNotFoundError("len_snapshots file do not exist!\n")
    
    n_rank = len( len_Y )
    
    modes_Y = copy.deepcopy( modes_temp )
    modes_Z = copy.deepcopy( modes_temp )
    
    for i in range( n_modes ):
        
        del modes_Y.modes[i].Phi
        del modes_Z.modes[i].Phi
        modes_Y.modes[i].Phi = []
        modes_Z.modes[i].Phi = []
        pos = 0
        
        for j in range( n_rank ):
            
            modes_Y.modes[i].Phi = modes_Y.modes[i].Phi + \
                              modes_temp.modes[i].Phi[pos:pos+len_Y[j]].tolist()
            pos = pos + len_Y[j]

            modes_Z.modes[i].Phi = modes_Z.modes[i].Phi + \
                              modes_temp.modes[i].Phi[pos:pos+len_Z[j]].tolist()
            pos = pos + len_Z[j] 
    
    # also distribute the reconstructed data
    
    pos = 0
    modes_Y.recons_std_dmd = []
    modes_Z.recons_std_dmd = []
    
    print(f"Got len_Y is {len_Y},len_Z is {len_Z}")
    
    for j in range( n_rank ):

        modes_Y.recons_std_dmd = modes_Y.recons_std_dmd + \
                            modes_temp.recons_std_dmd[pos:pos+len_Y[j]].tolist()
        pos = pos + len_Y[j]
        
        modes_Z.recons_std_dmd = modes_Z.recons_std_dmd + \
                            modes_temp.recons_std_dmd[pos:pos+len_Z[j]].tolist()
        pos = pos + len_Z[j]
    
        
# -- reconstruct data

    modes_Y.reconstruct( step )
    modes_Z.reconstruct( step )

# -- print selected modes properties:
    
    print(f"\nGot { n_modes } modes:")
    print(modes_Y.indxes)
    
    print(f"\nAlpha_pol of these modes:")
    print(modes_Y.alphas_pol)
    
    print(f"\n|Alpha_pol| of these modes:")
    print([abs(alpha_pol) for alpha_pol in modes_Y.alphas_pol])
    
    print(f"\nmu of these modes:")
    print(modes_Y.mus)
    
    print(f"\n|mu| of these modes:")
    print([abs(mu) for mu in modes_Y.mus])
    
    print(f"\nSt of these modes:")
    print(modes_Y.Sts)   
    
    print(f"\nIndexes of positive modes:")
    print( modes_Y.df_ind )
    
sys.stdout.flush()

# =============================================================================
# match mesh and reconstructed data( both std_dmd and spdmd ), each modes
# =============================================================================
    
# match mesh with data
with timer("\n - Match mesh "):

    print("modes_Y:\n")
    modes_Y.match_mesh( dfY, 'Y' )
    print(modes_Y.df_modes)

    print("\nmodes_Z:\n")
    modes_Z.match_mesh( dfZ, 'Z' )
    print(modes_Z.df_modes)

sys.stdout.flush()

# =============================================================================
# generate the grid that will be interpolated on
# =============================================================================

with timer("\n - Generate the interpolation grid "):
    
    xmin = modes_Y.df_modes['x'].min()  # Y plane's range in x is smaller
    xmax = 119.0
#    xmax = modes_Y.df_modes['x'].max()
    
    ymin = 0.040
#    ymin = 0.0
    ymax = 44.04
#    ymax = modes_Z.df_modes['y'].max()
    
    zmin = modes_Y.df_modes['z'].min()
    zmax = modes_Y.df_modes['z'].max()
    
    x = np.linspace( xmin, xmax, 451)
    y = np.linspace( ymin, ymax, 111)
    z = np.linspace( zmin, zmax, 121)
    
    modes_Y.grids_interp = np.meshgrid( x, z )
    modes_Z.grids_interp = np.meshgrid( x, y )
    
    

# =============================================================================
# interpolate and output figures
# =============================================================================

with timer("\n - Interpolate and output figures "):
    
    # enter ./modesfigures directory
    
    if os.path.exists("./modesfigures"): os.chdir("./modesfigures")
    
    else: os.mkdir("./modesfigures"); os.chdir("./modesfigures")
    
    # plot reconstructed flow field with selected spdmd modes
    
    header = "recons_00000"
    
    xx1, yy1, v1 = modes_Y.interp_recons( header,'Y' )
    xx2, yy2, v2 = modes_Z.interp_recons( header,'Z' )
        
    plot_combined_dmd_mode( (xx1,yy1), v1, 'y',
                            (xx2,yy2), v2, 'z',
                            filename='reconstructed_spdmd.png',
                            colorbar=True,
                            title='reconstructed spdmd modes')
    
    # plot reconstructed flow field with std dmd modes
    
    header = "recons_std_dmd"

    xx1, yy1, v1 = modes_Y.interp_recons( header,'Y' )
    xx2, yy2, v2 = modes_Z.interp_recons( header,'Z' )
    
    plot_combined_dmd_mode( (xx1,yy1), v1, 'y',
                            (xx2,yy2), v2, 'z',
                            filename='reconstructed_std_dmd.png',
                            colorbar=True,
                            title='reconstructed full modes')

    # indexes of positive modes

    indxes = np.array( modes_Y.df_ind['indxes'] )

    print(f"positive modes indexes: {indxes}")
      
    # plot mean mode
    
    print(f"The index of max alpha_pol is {indxes[0]}.\n")
    
    xx1, yy1, v1 = modes_Y.interp_mode( indxes[0], 'Y' )
    xx2, yy2, v2 = modes_Z.interp_mode( indxes[0], 'Z' )
    
    plot_combined_dmd_mode( (xx1,yy1), v1, 'y',
                            (xx2,yy2), v2, 'z',
                            filename='mean_mode.png',
                            colorbar=True,
                            title='mean mode')    
    
    # plot oscillating dmd modes
    
    phases = [cmath.rect(1.0, cmath.pi*0.25*i) for i in range(8)]
    
    Sts = np.array( modes_Y.df_ind['Sts'] )
    
    for n in range( 1, len(indxes) ):
        
        cmax = 0.0
        
        for i, phase in enumerate(phases):
            
            xx1, yy1, v1 = modes_Y.interp_mode( indxes[n], 'Y', phase=phase )
            xx2, yy2, v2 = modes_Z.interp_mode( indxes[n], 'Z', phase=phase ) 
        
            cmax_new = max(np.max(np.abs(v1)),np.max(np.abs(v2)))
        
            cmax = max(cmax, cmax_new)
        
        clevel = np.linspace( -cmax, cmax, 51 )
        
        for i, phase in enumerate(phases):
            
            xx1, yy1, v1 = modes_Y.interp_mode( indxes[n], 'Y', phase=phase )
            xx2, yy2, v2 = modes_Z.interp_mode( indxes[n], 'Z', phase=phase )         
    
            filename = f"mode_{n:02d}_{indxes[n]:05d}_{i}"
            
            title = f"St={Sts[n]:.3f} "+ f"phase = {i}/4"+r"$\pi$"
            
            plot_combined_dmd_mode( (xx1,yy1), v1, 'y',
                                    (xx2,yy2), v2, 'z',
                                    filename=filename+'.png',
                                    colorbar=True,
                                    clevel=clevel,
                                    title=title)
        
        # convert png image to gif
        
        convert_gif = f"convert -delay 100 -layers Optimize mode_{n:02d}*.png mode_{n:02d}.gif"
        
        exit_code = os.system( convert_gif )
        
        if exit_code == 0: print( f"Output mode_{n:02d} gif.\n")
        else:              print( f"Failed to output mode_{n:02d} gif.")
        

t_end = time.time()

print(f"\n - dmdpost totally took {t_end-t_0:8.2f}s.")


# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
sys.stdout.flush()       
    


#     header = "recons_std_dmd"
# #    header = "phi_00000"
    
#     v = np.array(modes_Y.df_modes[header]).real
    
#     x = np.unique(modes_Y.df_modes['x'])
#     z = np.unique(modes_Y.df_modes['z'])
    
#     n1 = len(x)
#     n2 = len(z)
    
#     xx = np.array(modes_Y.df_modes['x']).reshape(n2,n1)
#     zz = np.array(modes_Y.df_modes['z']).reshape(n2,n1)
    
#     v = v.reshape(n2,n1)
    
#     fig, ax = plt.subplots()
    
#     clevel = np.linspace( np.min(v),np.max(v),51)
    
#     cs =  ax.contourf(xx,zz,v,levels=clevel,cmap='bwr')
#     plt.colorbar(cs)
    
#     plt.show()

    
#     x = np.unique(modes_Z.df_modes['x'])
#     y = np.unique(modes_Z.df_modes['y'])
    
#     n1 = len(x)
#     n2 = len(y)
    
#     xx = np.array(modes_Z.df_modes['x']).reshape(n2,n1)
#     yy = np.array(modes_Z.df_modes['y']).reshape(n2,n1)
    
#     v = np.array(modes_Z.df_modes[header]).real
#     v = v.reshape(n2,n1)
    
#     fig, ax = plt.subplots()
    
#     clevel = np.linspace( np.min(v),np.max(v),51)
    
#     cs =  ax.contourf(xx,yy,v,levels=clevel,cmap='bwr')
#     plt.colorbar(cs)
    
#     plt.show()    