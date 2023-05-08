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

from   vista.paradmd     import ParaDmd


#snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots'
#
#snap_file = snap_dir+'/snap_test/snapshot_00600031/snapshot_W_002.bin'

snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots'

paradmd = ParaDmd( snap_dir )

paradmd.comm.barrier()

with timer('Read in '):
    
    paradmd.read_info_file()

    paradmd.read_struct_file()

    paradmd.assign_block()

    print('I am process %d out of %d, got the total block number %d.\
            I take care of blocks from %d(%d) .. %d(%d)'
            %(paradmd.rank, paradmd.n_procs, paradmd.n_bl,
            paradmd.i_start+1, paradmd.bl_num[paradmd.i_start], 
            paradmd.i_end, paradmd.bl_num[paradmd.i_end-1])) 

#    paradmd.comm.barrier()

    snap_files = None

    if paradmd.rank == 0:
        
        snap_files = get_filelist( snap_dir + '/snap_test' )
        
    snap_files = paradmd.comm.bcast( snap_files, root=0 )

    paradmd.select = 'p'

#    outputname = f'rank_{paradmd.rank}_list.dat'
#    
#    with open(outputname,'w') as f:
#        
#        for snap_file in snap_files:
#            
#            f.write(snap_file+'\n')
    for snap_file in snap_files:
        
        paradmd.para_read_data( snap_file )
        
    print(np.shape(paradmd.snapshots))

with timer('paradmd '):

    paradmd.do_paradmd()
