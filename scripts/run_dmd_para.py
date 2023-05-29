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



#snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots'
#
#snap_file = snap_dir+'/snap_test/snapshot_00600031/snapshot_W_002.bin'
#os.chdir( snap_dir )


snap_dir = os.getcwd()

paradmd = ParaDmd( snap_dir )


paradmd.comm.barrier()


with timer('Read in '):
    
    
    paradmd.read_info_file()

    paradmd.read_struct_file()

    paradmd.assign_block()

    print('I am process %d out of %d, got the total block number %d.\
            I take care of blocks from %d(%d) .. %d(%d)'
            %(paradmd.rank,      paradmd.n_procs,   paradmd.n_bl,
              paradmd.i_start+1, paradmd.bl_num[paradmd.i_start], 
              paradmd.i_end,     paradmd.bl_num[paradmd.i_end-1])) 

#    paradmd.comm.barrier()


    snap_files = None

    if paradmd.rank == 0:
        
        snap_files = get_filelist( snap_dir + '/snapshots' )
        
    snap_files = paradmd.comm.bcast( snap_files, root=0 )
    

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


