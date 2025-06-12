#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   slice_3dsnapshot_x.py
@Time    :   2024/11/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   slice 3d snapshot write out 2d slices
'''

import os
import sys
import time
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.mpi         import MPIenv
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.tools       import get_filelist
from   vista.snapshot    import Snapshot
from   vista.directories import Directories

# =============================================================================
# - build MPI communication environment

mpi = MPIenv()

# =============================================================================

casedir   = '/home/wencan/temp/smooth_mid'
slic_type = 'X'
locs      = [-14.98,-8.53,-7.33,1.89,8.0]
locs      = np.array(locs) * 5.2 + 50.4
bbox      = [-30.0,120.0,-1.0,31.0,-11.0,11.0]

# =============================================================================

dirs      = Directories( casedir )

# allocate variables that will be broadcasted to all workers

grid3d    = None
snapfiles = None

# root does the preparation

if mpi.is_root:
    
    # get all snapshot files
    
    snapfiles = get_filelist( dirs.snp_dir, 'snapshot.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot files")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    
grid3d    = mpi.comm.bcast( grid3d, root=0 )
snapfiles = mpi.comm.bcast( snapfiles, root=0 )

clock = timer(f"slice 3d snapshot in case {casedir}:")

# =============================================================================
# - action on one snapshot file

def slice_snapshot( snapfile ):
    
    # loop over all x locations
    
    for i, loc in enumerate(locs):
        
        outfile = os.path.dirname( snapfile) + f"/snapshot_{slic_type}_{i:03d}.bin"
        
        blocklist, _ = grid3d.select_sliced_blockgrids( slic_type, loc, bbox )
        
        snap = Snapshot( snapfile )
        snap.read_snapshot( blocklist )
        snap.grid3d = grid3d
        snap2d = snap.get_slice( slic_type, loc )
        snap2d.write_snapshot( outfile )
        
        print(f"Finish writing slice {snap2d.itstep:08d}_{i:03d}.")
        
        del snap, snap2d        
        

# =============================================================================

if mpi.size == 1:
    
    print("No worker available. Master should do all tasks.")
    
    for i, snapfile in enumerate(snapfiles):
        
        slice_snapshot( snapfile )
        clock.print_progress( i, len(snapfiles) )
    
else:
    if mpi.rank == 0:
        mpi.master_distribute( snapfiles )
    else:
        while True:
            task_index = mpi.worker_receive()
            
            if task_index is None: break
            else: 
                slice_snapshot( snapfiles[task_index] )
                clock.print_progress( task_index, len(snapfiles), rank=mpi.rank )

mpi.barrier()

if mpi.rank == 0:

    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush()     
