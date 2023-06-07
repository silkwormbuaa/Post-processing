# -*- coding: utf-8 -*-
'''
@File    :   vista_dmd_para.py
@Time    :   2023/05/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   parallel dmd script
'''

import os

import sys 

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

import numpy             as     np

import pandas            as     pd

from   mpi4py            import MPI

from   vista.timer       import timer

from   vista.tools       import get_filelist

from   vista.tools       import read_case_parameter

from   vista.paradmd     import ParaDmd

from   vista.snapshot    import Snapshot



#snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots'
#
#snap_file = snap_dir+'/snap_test/snapshot_00600031/snapshot_W_002.bin'
#os.chdir( snap_dir )


snap_dir = os.getcwd()

paradmd = ParaDmd( snap_dir )

# =============================================================================
# Get the snapshots info, struct files and case parameters, then broadcast.
# =============================================================================

snap_files = None

case_parameters = None

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

with timer('Read in all the snapshots'):
    
    paradmd.read_info_file()

    paradmd.read_struct_file()

    paradmd.assign_block()


    # verbose which processor take care of which blocks
    
    print(f"I am process {paradmd.rank:5d} out of {paradmd.n_procs:5d},",end='')
    print(f"got {paradmd.n_bl_local:5d} blocks out of {paradmd.n_bl:5d}",end='')
    print(f". I take care of blocks from {paradmd.i_start+1:5d}",end='')
    print(f"({paradmd.bl_num[paradmd.i_start]:5d}) to ",end='')
    print(f"{paradmd.i_end:5d}({paradmd.bl_num[paradmd.i_end-1]:5d}).")
    
    # select which parameter will be chosen to do DMD
    
    paradmd.select = 'p'
    paradmd.var_norms['p'] = float( case_parameters.get('p_ref') )
    
    for snap_file in snap_files:
        
        paradmd.para_read_data( snap_file )
    
    print(f"Finish read in snapshots.",end='')
    print(f"Snapshots shape of {paradmd.rank} is {np.shape(paradmd.snapshots)}")


# Specify the time interval of snapshots

paradmd.dt = float( case_parameters.get('dt_snap') )


with timer('paradmd '):

    paradmd.do_paradmd()
    
    # if spdmd is already done, just output selected modes

    if os.path.exists( paradmd.spdmd_result_file ):
        
        paradmd.read_spdmd_result()
        
        paradmd.para_write_modes()
        
    # if spdmd is yet to do
    
    else:
            
        if paradmd.rank == 0: paradmd.save_Pqs()
            


