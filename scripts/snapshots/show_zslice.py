#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   z_slice_DS.py
@Time    :   2024/02/20 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Visualize the DS contour from snapshot_Z_xxx.bin
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
from   vista.plot_style  import plot_slicez_stat
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
    
    if os.path.exists('./figures_u') == False: 
        os.system('mkdir ./figures_u')

    if os.path.exists('./figures_DS') == False: 
        os.system('mkdir ./figures_DS')
    
    if os.path.exists('./pkl_data') == False:
        os.system('mkdir ./pkl_data')

    snapfiles = get_filelist( snaps_dir, 'snapshot_Z' )

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

        u_slice = np.array( df_slice['u'] )
        grad_rho_slice = np.array( df_slice['grad_rho'] )
#        DS_slice = compute_DS( df_slice['grad_rho'] )

        x = np.linspace( -20, 12, 321 )

        y = np.linspace( 0.01, 8, 201 )

        xx,yy = np.meshgrid( x, y )

        u = griddata( (x_slice,y_slice), u_slice,
                      (xx,yy), method='linear')
        
        grad_rho = griddata( (x_slice,y_slice), grad_rho_slice,
                             (xx,yy), method='linear')
        
        DS = compute_DS( grad_rho, min=0.0, max=2.0 )
        
        with open(f'pkl_data/data_{snap.itstep:08d}.pkl', 'wb') as f:
            pickle.dump(xx, f)
            pickle.dump(yy, f)
            pickle.dump(u,  f)
            pickle.dump(DS, f)

        sep_line_file = f'separationlines_{snap.itstep:08d}.pkl'
        save_isolines( xx,yy,u, 0.0, sep_line_file )

        cbar = r'$u/u_{\infty}$'
        cbar_levels = np.linspace( -0.2, 1, 37)
        cbar_ticks  = np.linspace( -0.2, 1, 7)
        tag = f't = {snap.itime:6.2f} ms'
        plot_slicez_stat( xx,yy,u/507,
                          filename=f'u_{snap.itstep:08d}',
                          col_map='coolwarm',
                          cbar_label=cbar,
                          separation=sep_line_file,
                          sonic=False,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          tag=tag,
                          tag_loc=[-13,6],
                          x_lim=[-15,10],
                          y_lim=[0,8],
                          pure=False)
        
        cbar = r'$DS$'
        cbar_levels = np.linspace( 0.0, 0.8,33)
        
        plot_slicez_stat( xx,yy,DS,
                          filename=f'DS_{snap.itstep:08d}',
                          col_map='Greys_r',
                          cbar_label=cbar,
                          separation=sep_line_file,
                          sonic=False,
                          cbar_levels=cbar_levels,
                          tag=tag,
                          tag_loc=[-13,6],
                          x_lim=[-15,10],
                          y_lim=[0,8],
                          pure=False)
        
        
        os.system(f'mv u_{snap.itstep:08d}.png ./figures_u/')
        os.system(f'mv DS_{snap.itstep:08d}.png ./figures_DS/')
        os.system(f'mv {sep_line_file} ./pkl_data/')
        
        print(f"Processor {rank} finished {i+1}/{len(snapfiles)}.")
        sys.stdout.flush()

        