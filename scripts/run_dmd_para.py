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

source_dir = os.path.dirname( os.path.dirname( os.path.realpath(__file__) ))
sys.path.append( source_dir )

import numpy             as     np

import pandas            as     pd

from   mpi4py            import MPI

from   vista.timer       import timer

from   vista.tools       import get_filelist

from   vista.plot_style  import plot_eigens

from   vista.plot_style  import plot_amp_st

from   vista.paradmd     import ParaDmd

from   vista.snapshot    import Snapshot



#snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots'
#
#snap_file = snap_dir+'/snap_test/snapshot_00600031/snapshot_W_002.bin'
#os.chdir( snap_dir )


snap_dir = os.getcwd()

paradmd = ParaDmd( snap_dir )


# Get the snapshots info and struct files, then broadcast.

snap_files = None

with timer('Get snapshots file list and snapshots info'):
    
    if paradmd.rank == 0:
        
        snap_files = get_filelist( snap_dir + '/snapshots' )
        
        testfile = snap_files[0]
        
        snapshot_temp = Snapshot( testfile )
        
        snapshot_temp.get_snapshot_struct()
        
            
    snap_files = paradmd.comm.bcast( snap_files, root=0 )


paradmd.comm.barrier()


with timer('Read in '):
    
    
    paradmd.read_info_file()

    paradmd.read_struct_file()

    paradmd.assign_block()

    print(f"I am process {paradmd.rank:5d} out of {paradmd.n_procs:5d},",end='')
    print(f"got {paradmd.n_bl_local:5d} blocks out of {paradmd.n_bl:5d}",end='')
    print(f". I take care of blocks from {paradmd.i_start+1:5d}",end='')
    print(f"({paradmd.bl_num[paradmd.i_start]:5d}) to ",end='')
    print(f"{paradmd.i_end:5d}({paradmd.bl_num[paradmd.i_end-1]:5d}).")
    

#    paradmd.comm.barrier() 

    paradmd.select = 'p'
    paradmd.var_norms['p'] = 45447.289

    for snap_file in snap_files:
        
        paradmd.para_read_data( snap_file )
        
    print(np.shape(paradmd.snapshots))


# Specify the time interval of snapshots

paradmd.dt = 0.005

with timer('paradmd '):

    paradmd.do_paradmd()

    if paradmd.rank == 0:
        
        paradmd.save_Pqs()


