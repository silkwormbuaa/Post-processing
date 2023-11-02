#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   run_show_snap.py
@Time    :   2023/06/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import time

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

import numpy             as     np

import pandas            as     pd

from   vista.timer       import timer

from   vista.snapshot    import Snapshot

from   vista.paradmd     import ParaDmd

from   vista.tools       import get_filelist
from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_dmd_mode 

from   vista.colors      import colors    as col

from   scipy.interpolate import griddata

from   vista.log         import Logger
sys.stdout = Logger()

t_0 = time.time()
# =============================================================================
# read in one snapshot file to get grid vectors
# =============================================================================

#snap_dir = os.getcwd()

snap_dir = '/home/wencanwu/my_simulation/temp/220825_low/snapshots_W_test'
snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots/snapshot_test'
snap_dir = '/home/wencanwu/my_simulation/temp/220926_lowRe/snapshots/snapshot_test_W'
os.chdir(snap_dir)

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
    
    print(np.shape(GX))

    df = pd.DataFrame( GX, columns = GX_header )
    
    print(np.max( GX[:,0]),np.min(GX[:,0]))
    print(np.max( GX[:,1]),np.min(GX[:,1]))

paradmd = ParaDmd( snap_dir )
    
# =============================================================================
# Get the snapshots info, struct files and case parameters, then broadcast.
# =============================================================================

snap_files = None

case_parameters = None

print(col.bg.green,col.fg.red,"This is rank ",f"{paradmd.rank}",col.reset)

with timer('Get snapshots file list, snapshots info and case parameters'):
    
    if paradmd.rank == 0:
        
        snap_files = get_filelist( snap_dir + '/snapshots' )
        
        testfile = snap_files[0]
        
        snapshot_temp = Snapshot( testfile )
        
        snapshot_temp.get_snapshot_struct()
        
        case_parameters = read_case_parameter( 'case_parameters' )
        
            
    snap_files = paradmd.comm.bcast( snap_files, root=0 )
    
    case_parameters = paradmd.comm.bcast( case_parameters, root=0 )


paradmd.comm.barrier()

# =============================================================================
# Read in all the snapshots
# =============================================================================

with timer('\nRead in all the snapshots'):
    
    paradmd.read_info_file()

    paradmd.read_struct_file()

    paradmd.assign_block()


    # verbose which processor take care of which blocks
    
    print(f"I am process {paradmd.rank:5d} out of {paradmd.n_procs:5d},",end='')
    print(f"got {paradmd.n_bl_local:5d} blocks out of {paradmd.n_bl:5d}",end='')
    print(f". I take care of blocks from {paradmd.i_start+1:5d}",end='')
    print(f"({paradmd.bl_num[paradmd.i_start]:5d}) to ",end='')
    print(f"{paradmd.i_end:5d}({paradmd.bl_num[paradmd.i_end-1]:5d}).\n")
    
    # select which parameter will be chosen to do DMD
    
    paradmd.select = 'p'
    paradmd.var_norms['p'] = float( case_parameters.get('p_ref') )
    
    for snap_file in snap_files:
        
        paradmd.para_read_data( snap_file )
    
    print(f"Rank {paradmd.rank} finish read in snapshots ",end='')
    print(f"with shape of {np.shape(paradmd.snapshots)}")


# Specify the time interval of snapshots

paradmd.dt = float( case_parameters.get('dt_snap') )

snap1 = np.array(paradmd.snapshots[0])

snap2 = np.array(paradmd.snapshots[-1])

df['snap1'] = snap1

print(np.max(snap1),np.min(snap1))

print( len(snap1) )

df['snap2'] = snap2

x_1 = np.linspace(-50.0,100.0,101)
x_2 = np.linspace(-10.0,10.0,51)

meshgrid = np.meshgrid(x_1,x_2)

v =   griddata( (df['x'], df['z']),
                df['snap2'],
                (meshgrid[0],meshgrid[1]),
                method='linear' )

print(np.min(v),np.max(v))

plot_dmd_mode(meshgrid,v,
              filename='debug2.png',
              colorbar=True)

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()    




