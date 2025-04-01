#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   extract_pressure_at_pfmax.py
@Time    :   2025/04/01
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
point      = [0.0,0.12,0.01]
vars       = ['u','v','w','p','T']

# =============================================================================

dirs      = Directories( casefolder )

params    = Params( dirs.case_para_file )
x_pfmax   = params.x_pfmax*params.delta_0 + params.x_imp
point     = [x_pfmax,point[1],point[2]]

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
snap_steps = np.zeros( n_snaps,             dtype=int   )
snap_times = np.zeros( n_snaps,             dtype=float )
data       = np.zeros( (n_snaps,len(vars)), dtype=float )

# - extract the spanwise averaged pressure at pfmax location

def extract_pressure_at_pfmax( snapfile ):
    
    index = snapfiles.index( snapfile )
    bl_num, indx_probe = grd.find_probe_index( point )

    snap = Snapshot( snapfile )
    snap.grid3d=grd
    snap.read_snapshot( [bl_num], var_read=vars )

    data[index,:]     = snap.get_point_probed( point, vars )
    snap_times[index] = snap.itime
    snap_steps[index] = snap.itstep
    
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
data       = mpi.comm.reduce( data,       root=0, op=mpi.MPI.SUM )

# - write out the results

if mpi.is_root:
    
    os.chdir( create_folder( dirs.pp_snp_pfmax ) )

    with open('point_probe.dat','w') as f:
        f.write("itstep".rjust(10)+"itime".rjust(15))
        for var in vars:
            f.write( str(var).rjust(18) )
        f.write("\n")
                
        for i in range(n_snaps):
            f.write(f"{snap_steps[i]:10d}{snap_times[i]:15.3f}")
            f.write(f"{''.join(f'{x:18.4f}' for x in data[i,:])}\n")
        print(f"point_probe.dat is saved in {dirs.pp_snp_pfmax}.")
        
    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}.")
    sys.stdout.flush()
