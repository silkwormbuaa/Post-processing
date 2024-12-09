#!/usr/bin/env python3
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
import time
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.timer       import timer
from   vista.colors      import colors    as col
from   vista.params      import Params
from   vista.paradmd     import ParaDmd
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.directories import create_folder
#from   vista.log         import Logger
#sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

casepath = '/home/wencan/temp/231124'

dirs = Directories( casepath )

snap_dir = dirs.snp_dir

os.chdir( create_folder(dirs.pp_dmd) )

paradmd = ParaDmd( snap_dir )

# =============================================================================
# Get the snapshots info, struct files and case parameters, then broadcast.
# =============================================================================

snap_files = None

params = None

print(col.bg.green,col.fg.red,"This is rank ",f"{paradmd.rank}",col.reset)
sys.stdout.flush()

with timer('Get snapshots file list, snapshots info and case parameters'):
    
    if paradmd.rank == 0:
        
        snap_files = get_filelist( snap_dir, 'snapshot_X_004.bin' )
        
        testfile = snap_files[0]
        
        snapshot_temp = Snapshot( testfile )
        
        snapshot_temp.get_snapshot_struct()
        
        params = Params( dirs.case_para_file )
            
    snap_files = paradmd.comm.bcast( snap_files, root=0 )
    
    params = paradmd.comm.bcast( params, root=0 )

sys.stdout.flush()

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
    sys.stdout.flush()
    
    # select which parameter will be chosen to do DMD
    
    paradmd.select = ['u','v','w','p']
    paradmd.var_norms['p'] = params.p_ref
    paradmd.var_norms['u'] = params.u_ref
    paradmd.var_norms['v'] = params.u_ref
    paradmd.var_norms['w'] = params.u_ref
    
    for snap_file in snap_files:
        
        paradmd.para_read_data( snap_file )
    
    print(f"Rank {paradmd.rank} finish read in snapshots ",end='')
    print(f"with shape of {np.shape(paradmd.snapshots)}")
    sys.stdout.flush()

# Specify the time interval of snapshots

paradmd.dt = params.dt_snap

# =============================================================================
# do parallel dmd
# =============================================================================

with timer('paradmd '):

    paradmd.do_paradmd()
    
    # if spdmd is already done, just output selected modes

    if os.path.exists( paradmd.spdmd_result_file ):
        
        paradmd.read_spdmd_result()
        
        paradmd.para_write_modes()
        
        os.chdir( dirs.pp_dmd )
        paradmd.save_reconstruct()
                
    # if spdmd is yet to do
    
    else:
            
        if paradmd.rank == 0: paradmd.save_Pqs()

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()       
