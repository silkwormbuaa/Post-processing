#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   separation_pdf.py
@Time    :   2024/05/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Script of computing the possibility density function of separation.

             Need inca_grid.bin and snapshots/ in the same directory.
'''


import os
import sys
import pickle
import numpy             as     np
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.directories import create_folder
from   vista.timer       import timer

comm    = MPI.COMM_WORLD
rank    = comm.Get_rank()
n_procs = comm.Get_size()

# =============================================================================

casefolder = '/home/wencan/temp/250304'

# =============================================================================

dirs          = Directories( casefolder )

snapshotfiles = None

print(f"Rank {rank} is working in folder: {dirs.case_dir}.")
sys.stdout.flush()

# root processor get the filelist

if rank == 0:
    
    print(f"I am root, now at {dirs.case_dir}.")
    sys.stdout.flush()
    
    snapshotfiles = get_filelist( dirs.snp_dir, key='snapshot.bin')

# broadcast

snapshotfiles = comm.bcast( snapshotfiles, root=0)

# distribute the tasks

n_snaps = len( snapshotfiles )
i_start, i_end = distribute_mpi_work( n_snaps, n_procs, rank )
print(f"Rank {rank} is working on snapshots from {i_start} to {i_end}.")
sys.stdout.flush()

# using one snapshot as a container to store separation times

snap_container = Snapshot( snapshotfiles[i_start] )
snap_container.read_snapshot( var_read=['u'] )

for bl in snap_container.snap_data:
    
    bl.df['n_sep'] = np.zeros_like( bl.df['u'], dtype=int ) 

# loop over the snapshots to count the separation times

for i, snapshotfile in enumerate( snapshotfiles[i_start:i_end] ):
    
    with timer(f'load snapshot data ...{snapshotfile[-25:]}'):
        
        snap = Snapshot( snapshotfile )
        snap.read_snapshot( var_read=['u'] )
    sys.stdout.flush()
    
    for j, bl in enumerate(snap.snap_data):
        
        bl_container = snap_container.snap_data[j]
        
        # count the separation times
        
        bl_container.df['n_sep'] += (bl.df['u'] < 0.0)
    

with timer("communication"):
    
    for bl in snap_container.snap_data:
        
        bl.df['n_sep']   = comm.reduce( np.array(bl.df['n_sep']), op=MPI.SUM, root=0)
        bl.df['pdf_sep'] = bl.df['n_sep'] / len(snapshotfiles) 

comm.barrier()

if rank == 0:
    os.chdir( create_folder(dirs.pp_bubble) )
    with open('snapshot_container.pkl','wb') as f:
        pickle.dump( snap_container, f )

