#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   dmd_para.py
@Time    :   2025/06/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.params      import Params
from   vista.paradmd     import ParaDmd
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.directories import create_folder


def main():

    case_dir = '/home/wencan/temp/220927/'
    dirs     = Directories( case_dir )
    
    snap_dir = dirs.snp_dir
    
    paradmd  = ParaDmd( snap_dir )
    
    snapfiles_z = None
    snapfiles_y = None
    params      = None
    
    if paradmd.rank == 0:

        snapfiles_z = get_filelist( snap_dir, 'snapshot_Z_001' )
        snapfiles_y = get_filelist( snap_dir, 'snapshot_Y_003' )
        params      = Params( dirs.case_para_file )
            
    snap_files_z = paradmd.comm.bcast( snapfiles_z, root=0 )
    snap_files_y = paradmd.comm.bcast( snapfiles_y, root=0 )
    params       = paradmd.comm.bcast( params, root=0 )

    i_s, i_e = distribute_mpi_work( len(snap_files_z), paradmd.n_procs, paradmd.rank )
    
    data_lists = []
    
    for i in range( i_s, i_e ):

        data_vector = np.concatenate( [data_from_snapshots( snap_files_z[i],vars=['u'] ), 
                                       data_from_snapshots( snap_files_y[i],vars=['u'] )] )
        data_lists.append( data_vector )

    data = np.array( data_lists )
    
    print( data.shape )
    print( data.dtype )
    
    paradmd.non_blocking_data_exchange( data )
    
    paradmd.dt = params.dt_snap
    paradmd.do_paradmd()
    
    if os.path.exists( paradmd.spdmd_result_file ):
        
        paradmd.read_spdmd_result()
        paradmd.para_write_modes()
        
        os.chdir( dirs.pp_dmd )
    
    if paradmd.rank == 0:
        
        paradmd.save_Pqs()
    
        print("DMD analysis completed.")
    
        
def data_from_snapshots( snapfile, blocklist=None, vars=None ):
    
    snp = Snapshot( snapfile )
    snp.read_snapshot(var_read=vars, block_list=blocklist)
    snp.drop_ghost()
    snp.assemble_block()
    
    return snp.df[vars].values.ravel()
        
        

# =============================================================================
if __name__ == "__main__":

    main()

