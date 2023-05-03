# -*- coding: utf-8 -*-
'''
@File    :   vista_dmd_para.py
@Time    :   2023/05/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   parallel dmd script
'''

import os

import sys 

import numpy             as     np

import pandas            as     pd

from   mpi4py            import MPI

sys.path.append('..')

from   utils.timer       import timer

from   utils.tools       import get_filelist

from   vista.snapshot    import Snapshot


# ----------------------------------------------------------------------
# >>> Input area                                                
# ----------------------------------------------------------------------

snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots/snap_test'

# ==============================================================

# Setup the communicator for parallel computation

comm = MPI.COMM_WORLD

rank = comm.Get_rank()

n_procs = comm.Get_size()


# Reading the filelist from the root process

files = None

n_snap_total = None

if (rank == 0):
    
    files = get_filelist( snap_dir )
    
    n_snap_total = len( files )
    

# Broadcasting the total number of "snapshots n_snap_total" in root
# into variable "n_snap_total"

files = comm.bcast( files, root=0 )

n_snap_total = comm.bcast( n_snap_total, root=0 )

# Calculating the snapshots that I should take care: "n_snap"
# Calculating the start index "i_start" within 0..n_snap_total-1
# or i_start = None and n_snap = 0 if there is no snapshots

# The following algorithm divides snapshots in a dealing way
# (I mean the number of snapshots instead of the actual snapshots)

n_snap = n_snap_total // n_procs  # floor dividing, i.e. rounding off

left_procs = n_snap_total - n_procs*n_snap


if ( rank < left_procs ):
    
    n_snap = n_snap + 1
    
    i_start = 0 + rank*n_snap
    
    i_end = i_start + n_snap

    
elif ( n_snap > 0 ):
    
    i_start = 0 + rank*n_snap + left_procs
    
    i_end = i_start + n_snap


else:
    
    n_snap = 0
    
    i_start = None
    
    i_end = None
    


print(f'I am process {rank} out of {n_procs}, got the total snapshot number \
      {n_snap_total}. I take care of snapshots {i_start+1} .. {i_end}') 
