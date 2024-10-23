#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   run_dmd_para_combined.py
@Time    :   2023/08/25 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import time
import pickle
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.timer       import timer
from   vista.params      import Params
from   vista.paradmd     import ParaDmd
from   vista.snapshot    import Snapshot
from   vista.colors      import colors    as col
from   vista.tools       import get_filelist
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

snap_dir = os.getcwd()

os.chdir(snap_dir)

paradmd = ParaDmd( snap_dir )

paradmd2 = ParaDmd( snap_dir )

# =============================================================================
# Get the snapshots info, struct files and case parameters, then broadcast.
# =============================================================================

snap_filesY = None
snap_filesZ = None

case_parameters = None

print(col.bg.green,col.fg.red,"This is rank ",f"{paradmd.rank}",col.reset)

with timer('Get Y snapshots file list, snapshots info and case parameters'):
    
    if paradmd.rank == 0:
        
        snap_filesY = get_filelist( snap_dir + '/snapshots','snapshot_Y_002.bin')
        
        testfile = snap_filesY[0]
        
        snapshot_temp = Snapshot( testfile )
        
        snapshot_temp.get_snapshot_struct('snap_structY.csv','snap_infoY.dat')
        
        case_parameters = Params( 'case_parameters' )
        
            
    snap_filesY = paradmd.comm.bcast( snap_filesY, root=0 )
    
    case_parameters = paradmd.comm.bcast( case_parameters, root=0 )

    paradmd.comm.barrier()

with timer('Get Z snapshots file list, snapshots info and case parameters'):
    
    if paradmd2.rank == 0:
        
        snap_filesZ = get_filelist( snap_dir + '/snapshots','snapshot_Z_001.bin')
        
        testfile = snap_filesZ[0]
        
        snapshot_temp = Snapshot( testfile )
        
        snapshot_temp.get_snapshot_struct('snap_structZ.csv','snap_infoZ.dat')        
    
            
    snap_filesZ = paradmd2.comm.bcast( snap_filesZ, root=0 )

    paradmd2.comm.barrier()

    sys.stdout.flush()

# =============================================================================
# Read in all the snapshots
# =============================================================================

with timer('\nRead in all the snapshots Y'):
    
    paradmd.read_info_file('snap_infoY.dat')

    paradmd.read_struct_file('snap_structY.csv')

    paradmd.assign_block()


    # verbose which processor take care of which blocks
    
    print(f"I am process {paradmd.rank:5d} out of {paradmd.n_procs:5d},",end='')
    print(f"got {paradmd.n_bl_local:5d} blocks out of {paradmd.n_bl:5d}",end='')
    print(f". I take care of blocks from {paradmd.i_start+1:5d}",end='')
    print(f"({paradmd.bl_num[paradmd.i_start]:5d}) to ",end='')
    print(f"{paradmd.i_end:5d}({paradmd.bl_num[paradmd.i_end-1]:5d}).\n")
    
    # select which parameter will be chosen to do DMD
    
    paradmd.select = 'p'
    paradmd.var_norms['p'] = case_parameters.p_ref
    
    for snap_file in snap_filesY:
        
        paradmd.para_read_data( snap_file )
    
    print(f"Rank {paradmd.rank} finish read in snapshots Y ",end='')
    print(f"with shape of {np.shape(paradmd.snapshots)}")


with timer('\nRead in all the snapshots Z'):
    
    paradmd2.read_info_file('snap_infoZ.dat')

    paradmd2.read_struct_file('snap_structZ.csv')

    paradmd2.assign_block()


    # verbose which processor take care of which blocks
    
    print(f"I am process {paradmd2.rank:5d} out of {paradmd2.n_procs:5d},",end='')
    print(f"got {paradmd2.n_bl_local:5d} blocks out of {paradmd2.n_bl:5d}",end='')
    print(f". I take care of blocks from {paradmd2.i_start+1:5d}",end='')
    print(f"({paradmd2.bl_num[paradmd2.i_start]:5d}) to ",end='')
    print(f"{paradmd2.i_end:5d}({paradmd2.bl_num[paradmd2.i_end-1]:5d}).\n")
    
    # select which parameter will be chosen to do DMD
    
    paradmd2.select = 'p'
    paradmd2.var_norms['p'] = case_parameters.p_ref
    
    for snap_file in snap_filesZ:
        
        paradmd2.para_read_data( snap_file )
    
    print(f"Rank {paradmd2.rank} finish read in snapshots Z ",end='')
    print(f"with shape of {np.shape(paradmd2.snapshots)}")

sys.stdout.flush()

# Specify the time interval of snapshots

paradmd.dt = case_parameters.dt_snap

# =============================================================================
# keep the length of snapshots and pass paradmd2's snapshot Z to paradmd1
# =============================================================================
paradmd.N_t, paradmd.len_snap_local = np.shape( paradmd.snapshots )
paradmd2.N_t, paradmd2.len_snap_local = np.shape( paradmd2.snapshots )

len_Y = paradmd.comm.gather( paradmd.len_snap_local, root = 0 )
len_Z = paradmd2.comm.gather( paradmd2.len_snap_local, root = 0 )

print(f"length of Y snapshots on each rank is {len_Y}")
print(f"length of Z snapshots on each rank is {len_Z}")

if paradmd.rank == 0:
    
    with open( paradmd.snap_dir + '/len_snapshots.pkl','wb' ) as f:
        
        pickle.dump( len_Y, f )
        pickle.dump( len_Z, f )
        
        print("Write the length of snapshots on each rank to len_snapshots.pkl.\n")


paradmd.snapshots = np.concatenate( ( paradmd.snapshots,
                                      paradmd2.snapshots ), axis=1 )

print( f"Rank {paradmd.rank} combined snapshots with a resulted ",end='')
print( f"shape of {np.shape(paradmd.snapshots)}.\n" )
sys.stdout.flush()

# =============================================================================
# do parallel dmd
# =============================================================================

with timer(f'Rank {paradmd.rank} paradmd '):

    paradmd.do_paradmd()
    
    # if spdmd is already done, just output selected modes

    if os.path.exists( paradmd.spdmd_result_file ):
        
        paradmd.read_spdmd_result()
        
        paradmd.para_write_modes()
        
        paradmd.save_reconstruct()
                
    # if spdmd is yet to do
    
    else:
            
        if paradmd.rank == 0: paradmd.save_Pqs()

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()       





