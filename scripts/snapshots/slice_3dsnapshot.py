#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   slice_3dsnapshot.py
@Time    :   2023/08/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   slice_3dsnapshot -> interpolate_zslice -> zslice_from_pkl
'''


import os
import sys
import gc
import time
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )


# =============================================================================

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_procs = comm.Get_size()


slic_type   = 'Y'
loc         = 10.4
output_file = '/snapshot_Y_003.bin'

# --- Root gets all the files and broadcast to other processors

filelist = None

if rank == 0:
    filelist = get_filelist( '.', 'snapshot.bin' )

filelist = comm.bcast( filelist, root=0 )


# --- Distribute the tasks (evenly as much as possible)
    
n_snaps = len( filelist )

i_start, i_end = distribute_mpi_work(n_snaps, n_procs, rank)

filelist = filelist[i_start:i_end]

print(f"I am processor {rank}, I take below tasks:")

for file in filelist:
    
    print(file)

print("=========="); sys.stdout.flush()

comm.barrier()

# Root read in grid file and broadcast to other processors

grid3d = None
block_list = None

# ----- check if the grid file is available and read in grid then broadcast

if rank == 0:
        
    if not os.path.exists('inca_grid.bin'):
        raise FileNotFoundError("Please check inca_grid.bin!")
    
    else:
        
        grid3d = GridData('inca_grid.bin')
        grid3d.verbose = False
        grid3d.read_grid()
        block_list, indx_slic = grid3d.select_sliced_blockgrids(slic_type,loc)
        
    sys.stdout.flush()
    
grid3d = comm.bcast( grid3d, root=0 )
block_list = comm.bcast( block_list, root=0 )


# ----- Slice snapshot one by one

i = 0

slice_size = None

for snapshot_file in filelist:
    
    # check if the snapshot slice already exists
    # if exists and slice file is complete(size equals), do not do slice
    
    do_slice = True
    
    slicefile = os.path.dirname(snapshot_file)+output_file
    
    if i > 0:
    
        file_is_there = os.path.exists(slicefile)
        
        if file_is_there: 
            
            slice_size = os.stat( slicefile ).st_size
            
            if slice_size == size_std: do_slice = False

    if do_slice:
    
        with timer("reading one snapshot and writing out slice"):

            # read in snapshot
            
            snapshot3d = Snapshot( snapshot_file )
            
            snapshot3d.verbose = False
            
            snapshot3d.read_snapshot( block_list )
            
            # write the slice
        
            snapshot3d.grid3d = grid3d
            
            snapshot2d = snapshot3d.get_slice( slic_type, loc )
                    
#            print(f"writing {slicefile}")
            
            snapshot2d.write_snapshot( slicefile )
            
            print(f"finish writing slice {snapshot2d.itstep}.")
            
            del snapshot3d,snapshot2d
            gc.collect()
    
    if i == 0: size_std = os.stat( slicefile ).st_size
    
    i += 1
    
    print(f"Rank {rank } progress: {i/len(filelist)*100:10.2f}%.",end='')
    print(f"  Finish file {snapshot_file[-30:]}")
        

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()     


# get a slice from a known folder

""" 
snapdir = '/home/wencanwu/my_simulation/temp/220926_lowRe/snapshots/snapshot_00452401/snapshot_block'

os.chdir( snapdir )

snapshot3d = Snapshot( 'snapshot.bin' )

snapshot3d.verbose = True

with timer("\nread 3d snapshot and get snapshot struct"):
    
    snapshot3d.get_snapshot_struct()

# write snapshot

with timer("\nwrite 3d snapshot"):
    
    snapshot3d.write_snapshot("snapshot3d.bin")

# get slice of snapshot3d

with timer("\ndo slice"):
    
    snapshot3d.grid3d = GridData('inca_grid.bin')
    snapshot3d.grid3d.verbose = True
    snapshot3d.grid3d.read_grid()
        
    snapshot2d = snapshot3d.get_slice( 'Y', 0.0 )
    
    snapshot2d.verbose = True
    
# write snapshots
    
with timer("\nwrite 2d snapshot"):
    
    snapshot2d.write_snapshot("snapshot2d.bin") 
 """

"""
# just show a slice that already there; !! only applicable to uniform mesh.

#os.chdir("/home/wencanwu/my_simulation/temp/220926_lowRe/snapshots/snapshot_00452401/snapshot_block")

with timer("\nread in new snapshot and show"):
    
    snapshot2d_new = Snapshot("snapshot_Y_003.bin")
    
    snapshot2d_new.verbose = False
    
    snapshot2d_new.read_snapshot()
    
    snapshot2d_new.drop_ghost( buff=3 )
    
    snapshot2d_new.assemble_block()
    
    print(snapshot2d_new.df)
    
    x = np.array( snapshot2d_new.df['x'] )
    
    z = np.array( snapshot2d_new.df['z'] )
    
    p = np.array( snapshot2d_new.df['p'] )

    N_x = len(np.unique(x))
    N_z = len(np.unique(z))
    
    x = x.reshape(N_z,N_x)
    z = z.reshape(N_z,N_x)
    p = p.reshape(N_z,N_x)
    
    fig, ax = plt.subplots()
    contour = ax.pcolor(x,z,p)
    ax.set_title('pressure')
    plt.show()
"""



