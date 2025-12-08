#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   zslice_pkl.py
@Time    :   2024/02/28 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   slice_3dsnapshot -> interpolate_zslice -> zslice_from_pkl
'''

import os
import sys
import pickle
import numpy             as     np
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.plot_style  import plot_slicez_stat
from   vista.plane_analy import compute_DS
from   vista.timer       import timer

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_procs = comm.Get_size()

case_folder = '/home/wencan/temp/smooth_mid'
dirs        = Directories(case_folder)

pklfolder = case_folder + '/postprocess/temporary/pkl'
os.chdir( pklfolder )

# --- distribute tasks

pklfiles = None

if rank == 0:
    print(f"I am root, now at {pklfolder}.")
    pklfiles = get_filelist( pklfolder, 'data_' )
    create_folder( pklfolder + '/../figures_u/' )
    create_folder( pklfolder + '/../figures_DS/')

pklfiles = comm.bcast( pklfiles, root=0 )
n_pklfiles = len( pklfiles )

i_start, i_end = distribute_mpi_work(n_pklfiles, n_procs, rank)
pklfiles = pklfiles[i_start:i_end]

print(f"I am processor {rank}, I take {i_end-i_start} tasks from {i_start} to {i_end}.")
sys.stdout.flush()

# --- process tasks

with timer("Plotting"):
        
    for file in pklfiles:
        with open(file, 'rb') as f:
            snapstep = pickle.load(f)
            snaptime = pickle.load(f)
            xx       = pickle.load(f)
            yy       = pickle.load(f)
            u        = pickle.load(f)
            p        = pickle.load(f)
            rho      = pickle.load(f)
            grad_rho = pickle.load(f)
    
        sep_line_file = f'separationlines_{snapstep:08d}.pkl'
    
        cbar = r'$u/u_{\infty}$'
        cbar_levels = np.linspace( -0.2, 1, 37)
        cbar_ticks  = np.linspace( -0.2, 1, 7)
        tag = f''
    
    
        plot_slicez_stat( xx,yy,u/507,
                        filename=f'u_{snapstep:08d}',
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
        DS = compute_DS( grad_rho, min=0.0, max=2.0 )
    
        plot_slicez_stat( xx,yy,DS,
                    filename=f'DS_{snapstep:08d}',
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
    
    
        os.system(f'mv u_{snapstep:08d}.png ../figures_u/')
        os.system(f'mv DS_{snapstep:08d}.png ../figures_DS/')
        sys.stdout.flush()
