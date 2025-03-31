#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   extract_pressure_at_pfmax.py
@Time    :   2025/03/31 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Extract the spanwise avearaged pressure at the pfmax location.
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

casefolder = '/home/wencan/temp/smooth_mid/'

# =============================================================================

dirs      = Directories( casefolder )

params    = Params( dirs.case_para_file )
grd       = None
snapfiles = None

print(f"Rank {mpi.rank} is working on {dirs.case_dir}.")
sys.stdout.flush()

# - root processor read in grid

if mpi.is_root:
    
    print(f"I am root, now at {dirs.case_dir}.")
    sys.stdout.flush()
    
    snapfiles = get_filelist( dirs.snp_dir, key='snapshot.bin' )
    
    with timer('load grid data'):
        grd = GridData( dirs.grid )
        grd.read_grid()
    sys.stdout.flush()
    
# - broadcast grid data to all processors

mpi.comm.barrier()
grd       = mpi.comm.bcast( grd,       root=0 )
snapfiles = mpi.comm.bcast( snapfiles, root=0 )

n_snaps    = len(snapfiles)
p_at_pfmax = np.zeros( n_snaps, dtype=float )
snap_steps = np.zeros( n_snaps, dtype=int )
snap_times = np.zeros( n_snaps, dtype=float )

# - extract the spanwise averaged pressure at pfmax location

def extract_pressure_at_pfmax( snapfile ):
    
    index = snapfiles.index( snapfile )
    loc   = [params.x_pfmax*params.delta_0 + params.x_imp, 0.001]
    blocklist, indx_probe = grd.select_probed_blockgrids( 'Z', loc )

    snap = Snapshot( snapfile )
    snap.read_snapshot( blocklist, var_read=['p'] )

    df_snap   = snap.get_probed_df( blocklist, grd, indx_probe, 'Z' )
    p_spanave = df_snap['p'].mean()
    
    snap_steps[index] = snap.itstep
    snap_times[index] = snap.itime
    p_at_pfmax[index] = p_spanave
    

# - distribute works

clock = timer("get pressure at pfmax from snapshots:")

if mpi.size == 1:
    print("No workers available. Master should do all tasks.")
    
    for i, snapfile in enumerate(snapfiles):
        extract_pressure_at_pfmax( snapfile )
        clock.print_progress( i, len(snapfiles), rank=mpi.rank )

else:
    if mpi.is_root:
        mpi.master_distribute( snapfiles )
    else:
        while True:
            task_index = mpi.worker_receive()
            
            if task_index is None: break
            else:
                extract_pressure_at_pfmax( snapfiles[task_index] )
                clock.print_progress( task_index, len(snapfiles), rank=mpi.rank )
                
# - gather results

mpi.barrier()

snap_steps = mpi.comm.reduce( snap_steps, root=0, op=mpi.MPI.SUM )
snap_times = mpi.comm.reduce( snap_times, root=0, op=mpi.MPI.SUM )
p_at_pfmax = mpi.comm.reduce( p_at_pfmax, root=0, op=mpi.MPI.SUM )

# - write out the results

if mpi.is_root:
    
    os.chdir( create_folder( dirs.pp_snp_pfmax ))

    with open('pressure_at_pfmax.dat','w') as f:
        f.write("itstep".rjust(10)+"itime".rjust(15))
        f.write("p_spanave_pfmax".rjust(20)+"\n")
                
        for i in range(n_snaps):
            f.write(f"{snap_steps[i]:10d}{snap_times[i]:15.3f}")
            f.write(f"{p_at_pfmax[i]:20.4f}\n")
        print(f"pressure_at_pfmax.dat is saved in {dirs.pp_snp_pfmax}.")
        
    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}.")
    sys.stdout.flush()
