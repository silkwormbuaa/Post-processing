#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   test_mpi.py
@Time    :   2024/09/12 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData

comm    = MPI.COMM_WORLD
rank    = comm.Get_rank()
n_procs = comm.Get_size()

a = GridData()

print(f"Hello World, I am rank {rank} out of {n_procs}.\n")
