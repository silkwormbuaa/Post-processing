#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   z_slice_DS.py
@Time    :   2024/02/20 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Interpolate z slice type snapshot and dump the data into pkl files
             for faster plotting.
'''

import os
import sys
import pickle
import numpy             as     np
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.plane_analy import shift_coordinates
from   vista.plane_analy import compute_DS
from   vista.plane_analy import save_isolines
from   scipy.interpolate import griddata
from   vista.timer       import timer

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_procs = comm.Get_size()

#snaps_dir = '/home/wencanwu/my_simulation/temp/220927_lowRe/snapshots/video_test/snapshots'

# -- get snapshots list

snaps_dir = os.getcwd()

os.chdir( snaps_dir )

snapfiles = None

if rank == 0:
    print(f"I am root, now at {snaps_dir}.")
    if os.path.exists('./pkl_data') == False:
        os.system('mkdir ./pkl_data')
    snapfiles = get_filelist( snaps_dir, 'snapshot_Z' )

# wait for the root to create the directory
comm.barrier() 
os.chdir('./pkl_data')
snapfiles = comm.bcast( snapfiles, root=0 )

# --- Distribute the tasks (evenly as much as possible)
    
n_snaps = len( snapfiles )
i_start, i_end = distribute_mpi_work(n_snaps, n_procs, rank)
snapfiles = snapfiles[i_start:i_end]

print(f"I am processor {rank}, I take below tasks:")

for file in snapfiles:
    print(file)

print("=========="); sys.stdout.flush()

comm.barrier()

#for file in snapfiles: print( file )
with timer('show snapshots'):
    for i, snapfile in enumerate(snapfiles):
        
        snap = Snapshot( snapfile )

        snap.verbose = False

        snap.read_snapshot()
        snap.compute_gradients()
        snap.drop_ghost( buff=3 )
        snap.assemble_block()

        delta = 5.2
        ximp  = 50.4

        df_slice = shift_coordinates( snap.df, delta, 0., 0., ximp )

        x_slice = np.array( df_slice['xs'] )
        y_slice = np.array( df_slice['y_scale'] )

        u_slice        = np.array( df_slice['u'] )
        p_slice        = np.array( df_slice['p'] )
        rho_slice      = np.array( df_slice['rho'] )
        grad_rho_slice = np.array( df_slice['grad_rho'] )
#        DS_slice = compute_DS( df_slice['grad_rho'] )

        x = np.linspace( -20, 12, 321 )
        y = np.linspace( 0.01, 8, 201 )
        xx,yy = np.meshgrid( x, y )

        u = griddata( (x_slice,y_slice), u_slice, (xx,yy), method='linear')
        p = griddata( (x_slice,y_slice), p_slice, (xx,yy), method='linear')
        rho = griddata( (x_slice,y_slice), rho_slice, (xx,yy), method='linear')
        grad_rho = griddata( (x_slice,y_slice), grad_rho_slice, (xx,yy), method='linear')
        
        DS = compute_DS( grad_rho, min=0.0, max=2.0 )
        
        with open(f'data_{snap.itstep:08d}.pkl', 'wb') as f:
            pickle.dump(snap.itstep, f)
            pickle.dump(snap.itime, f)
            pickle.dump(xx, f)
            pickle.dump(yy, f)
            pickle.dump(u,  f)
            pickle.dump(p,  f)
            pickle.dump(rho,f)
            pickle.dump(grad_rho, f)

        sep_line_file = f'separationlines_{snap.itstep:08d}.pkl'
        save_isolines( xx,yy,u, 0.0, sep_line_file )
        
        print(f"Processor {rank} finished {i+1}/{len(snapfiles)}.")
        sys.stdout.flush()
