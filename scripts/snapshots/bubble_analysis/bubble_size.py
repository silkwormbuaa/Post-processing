#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size.py
@Time    :   2024/04/22 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   compute the instanteneous separation bubble size considering the cutcells.
             Note: case/supplements should be ready. 
'''

import os
import sys
import time
import pickle
import numpy             as     np
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.mpi         import MPIenv
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.directories import create_folder

# - build MPI communication environment

mpi = MPIenv()

# =============================================================================

casefolder  = '/home/wencan/temp/231124/'

# =============================================================================

dirs = Directories( casefolder )

params        = Params( dirs.case_para_file )
roughwall     = params.roughwall
grd           = None
wd_snap       = None
snapshotfiles = None
cc_df         = None

print(f"Rank {mpi.rank} is working on {dirs.case_dir}.")
sys.stdout.flush()

# root processor read in grid, wall distance snapshot and cutcell info

if mpi.is_root:
    
    print(f"I am root, now at {dirs.case_dir}.")
    sys.stdout.flush()

    if roughwall:
        wd_snap_file = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]
        
    snapshotfiles = get_filelist( dirs.snp_dir, key='snapshot.bin')

    with timer('load grid data'):
        grd = GridData( dirs.grid )
        grd.read_grid()
        grd.cell_volume()
    sys.stdout.flush()

    if roughwall:
        
        with timer('load wall distance snapshot data'):
            wd_snap = Snapshot( wd_snap_file )
            wd_snap.read_snapshot( var_read=['wd'] )
        sys.stdout.flush()
        
        with timer("read cutcell info"):
            
            cc_df = pd.read_csv( dirs.cc_setup, delimiter=r'\s+')
            cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                        inplace=True)
        sys.stdout.flush()
    
    else: cc_df = None


# broadcast information to all the processors

mpi.comm.barrier()
grd           = mpi.comm.bcast( grd,           root=0 )
snapshotfiles = mpi.comm.bcast( snapshotfiles, root=0 )

if roughwall:
    wd_snap   = mpi.comm.bcast( wd_snap,       root=0 )
    cc_df     = mpi.comm.bcast( cc_df,         root=0 )

# - prepare the variables to store the bubble size

n_snaps           = len( snapshotfiles )
snapshot_step     = np.zeros(n_snaps, dtype=int)
snapshot_time     = np.zeros(n_snaps, dtype=float)
bubble_volume     = np.zeros(n_snaps, dtype=float)
bubble_volume_thr = np.zeros(n_snaps, dtype=float) 

# - use one snapshot as a container to store separation times

snap_container = Snapshot( snapshotfiles[0] )
snap_container.read_snapshot( var_read=['u'] )
for bl in snap_container.snap_data:
    bl.df['n_sep'] = np.zeros_like( bl.df['u'], dtype=int )

# - the work that workers do

def count_bubble_size( snapshotfile ):

    # load snapshot data
    
    index = snapshotfiles.index( snapshotfile )
    snap  = Snapshot( snapshotfile )
    snap.read_snapshot( var_read=['u'] )
    if roughwall:
        snap.assign_wall_dist( wd_snap )
    
    sys.stdout.flush()  
    
    # compute bubble volume  

    vol1, vol2 = snap.compute_bubble_volume( grd, cc_df=cc_df, roughwall=roughwall, y_threshold=0.0 )
    print(f"snapshot ...{snapshotfile[-25:]} bubble size:{vol1:15.4f} ({vol2:15.4f})")

    snapshot_step[index]     = snap.itstep
    snapshot_time[index]     = snap.itime
    bubble_volume[index]     = vol1
    bubble_volume_thr[index] = vol2
    sys.stdout.flush()
    
    # count the separation times
    
    for j, bl in enumerate(snap.snap_data):
        
        bl_container = snap_container.snap_data[j]
        bl_container.df['n_sep'] += (bl.df['u'] < 0.0)


# - distribute works

clock = timer("count bubble size from snapshots:")

if mpi.size == 1:
    print("No workers available. Master should do all tasks.")
    
    for i, snapfile in enumerate(snapshotfiles):
        count_bubble_size( snapfile )
        clock.print_progress( i, len(snapshotfiles), rank=mpi.rank )
        
else:
    if mpi.rank == 0:
        mpi.master_distribute( snapshotfiles )
    else:
        while True:
            task_index = mpi.worker_receive()
            
            if task_index is None: break
            else: 
                count_bubble_size( snapshotfiles[task_index] )
                clock.print_progress( task_index, len(snapshotfiles), rank=mpi.rank )

# wait for all the processors to finish computing bubble size

mpi.barrier()

snapshot_step     = mpi.comm.reduce( snapshot_step,     root=0, op=mpi.MPI.SUM )
snapshot_time     = mpi.comm.reduce( snapshot_time,     root=0, op=mpi.MPI.SUM )
bubble_volume     = mpi.comm.reduce( bubble_volume,     root=0, op=mpi.MPI.SUM )
bubble_volume_thr = mpi.comm.reduce( bubble_volume_thr, root=0, op=mpi.MPI.SUM )

for bl in snap_container.snap_data:
    bl.df['n_sep']   = mpi.comm.reduce( np.array(bl.df['n_sep']), op=mpi.MPI.SUM, root=0)
    bl.df['pdf_sep'] = bl.df['n_sep'] / len(snapshotfiles)

# root processor write the bubble size to file

if mpi.is_root:
    
    os.chdir( create_folder(dirs.pp_bubble) )
    
    with open('bubble_size.dat','w') as f:
        f.write("itstep".rjust(10)+"itime".rjust(15))
        f.write("bubble_volume".rjust(15)+"bubble_volume_thr".rjust(20)+"\n")
                
        for i in range(n_snaps):
            f.write(f"{snapshot_step[i]:10d}{snapshot_time[i]:15.3f}")
            f.write(f"{bubble_volume[i]:15.4f}{bubble_volume_thr[i]:20.4f}\n")

        print(f"bubble size data is saved in {dirs.pp_bubble}.")
        
    with open('snapshot_container.pkl','wb') as f:
        pickle.dump( snap_container, f )
        
        print(f"separation p.d.f is saved in {dirs.pp_bubble}/snapshot_container.pkl.")
        
    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush()