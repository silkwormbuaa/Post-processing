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

    case_dir = '/home/wencan/temp/231124/'
    vars     = ['u','v','w','p']
    
    dirs     = Directories( case_dir )
    
    snap_dir = dirs.snp_dir
    
    paradmd  = ParaDmd( snap_dir )
    
    snapfiles_x = None
    snapfiles_z = None
    snapfiles_y = None
    params      = None
    
    if paradmd.rank == 0:
        
        snapfiles_x = get_filelist( snap_dir, 'snapshot_X_004.bin' )
        snapfiles_y = get_filelist( snap_dir, 'snapshot_Y_000.bin' )
        snapfiles_z = get_filelist( snap_dir, 'snapshot_Z_001.bin' )
        params      = Params( dirs.case_para_file )

    snapfiles_x = paradmd.comm.bcast( snapfiles_x, root=0 )            
    snapfiles_y = paradmd.comm.bcast( snapfiles_y, root=0 )
    snapfiles_z = paradmd.comm.bcast( snapfiles_z, root=0 )
    params      = paradmd.comm.bcast( params,      root=0 )

    i_s, i_e = distribute_mpi_work( len(snapfiles_z), paradmd.n_procs, paradmd.rank )
    
    data_lists = []
    
    for i in range( i_s, i_e ):

        data_vector = np.concatenate( [data_from_snapshots( snapfiles_x[i],vars=vars ),
                                       data_from_snapshots( snapfiles_y[i],vars=vars ), 
                                       data_from_snapshots( snapfiles_z[i],vars=vars )] )
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
        paradmd.save_reconstruct()
        
    else:
        
        if paradmd.rank == 0:
            paradmd.save_Pqs()
            print("DMD analysis completed.")


def data_from_snapshots( snapfile, blocklist=None, vars=None ):
    
    snap = Snapshot( snapfile )
    snap.read_snapshot(var_read=vars, block_list=blocklist)
    snap.drop_ghost(block_list=blocklist)
    
    data_snap = []
    for num in snap.bl_nums_clean:
        
        block = snap.snap_cleandata[snap.bl_nums_clean.index(num)]
        
        for var in vars:
            if var in ['u','v','w']:
                block.df[var] /= 507.0
            elif var == 'p':
                block.df[var] /= 45447.289
                
        data_snap.append( block.df[vars].values.T.ravel() )
    
    return np.concatenate([data for data in data_snap])
        
        

# =============================================================================
if __name__ == "__main__":

    main()

