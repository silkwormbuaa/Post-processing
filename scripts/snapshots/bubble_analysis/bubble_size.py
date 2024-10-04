#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size.py
@Time    :   2024/04/22 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   compute the instanteneous separation bubble size considering the cutcells
             
             Note:
             
             for roughwall case: wall_dist folder, cutcells_setup.dat, 
                                 inca_grid.bin, snapshots folder  are need.
             for smoothwall case: inca_grid.bin, snapshots folder are need.
             
             update 'roughwall' to True or False.
             
             change the length of spename in grid.py.
'''

import os
import sys
import numpy             as     np
import pandas            as     pd
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.log         import Logger
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.timer       import timer
#sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

roughwall = False
dirs  = '/path/to/the/case'

# =============================================================================

comm    = MPI.COMM_WORLD
rank    = comm.Get_rank()
n_procs = comm.Get_size()

grd           = None
wd_snap       = None
snapshotfiles = None
cc_df         = None

print(f"Rank {rank} is working on {dirs}.")
sys.stdout.flush()


# root processor read in grid, wall distance snapshot and cutcell info

if rank == 0:
    
    print(f"I am root, now at {dirs}.")
    sys.stdout.flush()

    grid_file = dirs + '/results/inca_grid.bin'
    ccfile    = dirs + '/supplements/cutcells_setup.dat'

    if roughwall:
        wd_snap_file  = get_filelist( dirs + '/wall_dist', key='snapshot.bin')[0]
        
    snapshotfiles = get_filelist( dirs + '/snapshots', key='snapshot.bin')

    with timer('load grid data'):
        grd = GridData( grid_file )
        grd.read_grid()
        grd.cell_volume()
    sys.stdout.flush()

    if roughwall:
        
        with timer('load wall distance snapshot data'):
            wd_snap = Snapshot( wd_snap_file )
            wd_snap.read_snapshot( var_read=['wd'] )
        sys.stdout.flush()
        
        with timer("read cutcell info"):
            
            cc_df = pd.read_csv( ccfile, delimiter=r'\s+')
            cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                        inplace=True)
        sys.stdout.flush()
    
    else: cc_df = None


# broadcast information to all the processors

comm.barrier()
grd           = comm.bcast( grd, root=0 )
snapshotfiles = comm.bcast( snapshotfiles, root=0 )

if roughwall:
    wd_snap   = comm.bcast( wd_snap, root=0 )
    cc_df     = comm.bcast( cc_df, root=0 )


# distribute works

n_snaps = len( snapshotfiles )
i_start, i_end = distribute_mpi_work( n_snaps, n_procs, rank )

snapshot_index = np.zeros(n_snaps, dtype=int)
snapshot_time  = np.zeros(n_snaps, dtype=float)
bubble_volume  = np.zeros(n_snaps, dtype=float)


# compute bubble size

for i,snapshotfile in enumerate(snapshotfiles[i_start:i_end]):
    
    with timer(f'load snapshot data ...{snapshotfile[-15:]}'):
        
        snap = Snapshot( snapshotfile )
        snap.read_snapshot( var_read=['u'] )
        if roughwall:
            snap.assign_wall_dist( wd_snap )
    sys.stdout.flush()
    
    with timer('compute bubble size'):
        
        vol = snap.compute_bubble_volume( grd, cc_df=cc_df, roughwall=roughwall )
        
        print(f"snapshot ...{snapshotfile[-25:]} bubble size:{vol:15.4f}")
    
    snapshot_index[i_start+i] = snap.itstep
    snapshot_time[i_start+i]  = snap.itime
    bubble_volume[i_start+i]  = vol
    sys.stdout.flush()
    
    
# wait for all the processors to finish computing bubble size

snapshot_index = comm.reduce( snapshot_index, root=0, op=MPI.SUM )
snapshot_time  = comm.reduce( snapshot_time,  root=0, op=MPI.SUM )
bubble_volume  = comm.reduce( bubble_volume,  root=0, op=MPI.SUM )


# root processor write the bubble size to file

if rank == 0:
    
    with open('bubble_size.dat','w') as f:
        
        f.write("itstep itime bubble_volume\n")
        for i in range(n_snaps):
            f.write(f"{snapshot_index[i]:10d}{snapshot_time[i]:15.3f}{bubble_volume[i]:15.4f}\n")

