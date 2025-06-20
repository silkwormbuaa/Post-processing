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

from   vista.grid        import GridData
from   vista.params      import Params
from   vista.paradmd     import ParaDmd
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.directories import create_folder


def main():

    case_dir = '/home/wencan/temp/smooth_mid/'
    vars     = ['u','v','w','p']
    
    dirs     = Directories( case_dir )
    
    snap_dir = dirs.snp_dir
    
    paradmd  = ParaDmd( snap_dir )
    
    snapfiles_x = snapfiles_y = snapfiles_z = \
    stat_x      = stat_y      = stat_z      = \
    params      = None
    
    if paradmd.rank == 0:
        
        snapfiles_x = get_filelist( snap_dir, 'snapshot_X_004.bin' )
        snapfiles_y = get_filelist( snap_dir, 'snapshot_Y_000.bin' )
        snapfiles_z = get_filelist( snap_dir, 'snapshot_Z_001.bin' )

        grid        = GridData( dirs.grid )
        grid.read_grid()

        blocklist,_ = grid.select_sliced_blockgrids( 'X', 8.0*5.2+50.4 )
        stat_x      = read_statistic_slice( dirs.sup_dir + '/stat_xslice_004.bin', 
                                            blocklist, grid, vars=vars )
        blocklist,_ = grid.select_sliced_blockgrids( 'Y', 0.0001 )
        stat_y      = read_statistic_slice( dirs.sup_dir + '/stat_yslice_000.bin', 
                                            blocklist, grid, vars=vars )
        blocklist,_ = grid.select_sliced_blockgrids( 'Z', 0.0001 )
        stat_z      = read_statistic_slice( dirs.sup_dir + '/stat_zslice.bin',     
                                            blocklist, grid, vars=vars )
        
        params      = Params( dirs.case_para_file )

    snapfiles_x = paradmd.comm.bcast( snapfiles_x, root=0 )            
    snapfiles_y = paradmd.comm.bcast( snapfiles_y, root=0 )
    snapfiles_z = paradmd.comm.bcast( snapfiles_z, root=0 )
    stat_x      = paradmd.comm.bcast( stat_x,      root=0 )
    stat_y      = paradmd.comm.bcast( stat_y,      root=0 )
    stat_z      = paradmd.comm.bcast( stat_z,      root=0 )
    params      = paradmd.comm.bcast( params,      root=0 )

    i_s, i_e = distribute_mpi_work( len(snapfiles_z), paradmd.n_procs, paradmd.rank )
    
    data_lists = []
    
    for i in range( i_s, i_e ):

        data_vector = np.concatenate( [data_from_snapshots( snapfiles_x[i],stat_x,vars=vars ),
                                       data_from_snapshots( snapfiles_y[i],stat_y,vars=vars ), 
                                       data_from_snapshots( snapfiles_z[i],stat_z,vars=vars )] )
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


def data_from_snapshots( snapfile, stat:StatisticData, blocklist=None, vars=None ):
    
    snap = Snapshot( snapfile )
    snap.read_snapshot(var_read=vars, block_list=blocklist)
    snap.drop_ghost(block_list=blocklist)
    
    data_snap = []
    for num in snap.bl_nums_clean:
        
        block = snap.snap_cleandata[snap.bl_nums_clean.index(num)]
        
        statbl = stat.bl_clean[stat.bl_nums_clean.index(num)] if stat.bl_clean else None
        
        for var in vars:
            if var == 'u':
                block.df[var] = (block.df[var] - np.array(statbl.df[var],dtype=np.float32))/507.0
            elif var in ['v', 'w']:
                block.df[var] = (block.df[var] - np.array(statbl.df[var],dtype=np.float32))/152.10
            elif var == 'p':
                block.df[var] = (block.df[var] - np.array(statbl.df[var],dtype=np.float32))/45447.289
                
        data_snap.append( block.df[vars].values.T.ravel() )
    
    return np.concatenate([data for data in data_snap])


def read_statistic_slice( stat_file, blocklist, grid, vars=None):
    
    stat = StatisticData( stat_file )
    stat.grid3d = grid
    stat.read_statistic(block_list=blocklist, vars_in=vars)
    stat.drop_ghost(block_list=blocklist)
    
    return stat
        

# =============================================================================
if __name__ == "__main__":

    main()

