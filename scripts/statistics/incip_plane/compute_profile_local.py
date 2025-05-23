#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   compute_profile_local.py
@Time    :   2025/05/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Do spanwise periodic average first, then compute the profile at a local point.
'''


import os
import sys
import numpy             as     np
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.params      import Params
from   vista.tools       import get_filelist
from   vista.directories import create_folder


def main():
   
    case_dirs = ['smooth_adiabatic','220927',
                 'smooth_mid','231124','241030']
    
    for case in case_dirs:
        
        compute_profile_local( '/home/wencan/temp/' + case, loc='upstream' )
        compute_profile_local( '/home/wencan/temp/' + case, loc='incip' )
    
    
def compute_profile_local( case_folder, loc='upstream' ): 
   
    dirs         = Directories( case_folder )
    params       = Params( dirs.case_para_file )

    grid         = GridData( dirs.grid )
    grid.read_grid()
    grid.cell_volume()
    
    vars_read   = ['u','v','w','p','pp','rho','mu','T','uu','vv','ww','uv']
    vars_output = ['u','v','w','p','pp','rho','mu','T','uu','vv','ww','uv']

    if   loc == 'upstream': 
        loc_x = -53.6
        os.chdir( create_folder(dirs.pp_profile_up) )
    elif loc == 'incip':    
        loc_x = params.x_incip
        os.chdir( create_folder(dirs.pp_profile_incip) )
    
    z_valley  = 0.5*params.D
    H         = params.H
    roughwall = params.roughwall
    
    if roughwall:
        cc_df = pd.read_csv( dirs.cc_setup, delimiter = r'\s+' )
        cc_df.drop( columns=['fax0','faz0','fax1','faz1', 'processor']
                            , inplace=True )
    
    xz_ridge  = np.array([loc_x, 0.001])
    xz_valley = np.array([loc_x, z_valley+0.001])
    
    bbox_ridge  = [loc_x-0.1,loc_x+0.1, 0.0,100, -20, 20]
    bbox_valley = [loc_x-0.1,loc_x+0.1, -H ,   100, -20, 20]
    
# =============================================================================

    def process_stat(xz, bbox, outfile):

        blocklist_full = grid.select_blockgrids(bbox, mode='overlap')
        print( blocklist_full )
        
        
        blocklist,_    = grid.select_probed_blockgrids('Y',xz, bbox=bbox, bbox_mode='overlap')
        
        stat_file      = dirs.statistics
        stat           = StatisticData(stat_file)
        stat.grid3d    = grid
        stat.read_statistic(blocklist_full,vars_in=vars_read)
        
        if roughwall:
            wdfile  = get_filelist(dirs.wall_dist,key='snapshot.bin')[0]
            wd_snap = Snapshot( wdfile )
            wd_snap.read_snapshot( blocklist, var_read=['wd'] )
            for num in blocklist:
                temp_df = cc_df[cc_df['block_number'] == num]
                wall_dist = np.array( wd_snap.snap_data[num-1].df['wd'] )
                grid.g[num-1].assign_vol_fra( df=temp_df, wall_dist=wall_dist )
        else:
            for num in blocklist:
                grid.g[num-1].assign_vol_fra()
        
        stat.match_grid( blocklist_full, grid, add_to_df=True )
                
        # periodic averaging
        stat.spanwise_periodic_average( blocklist_full, vars_output, params.D )
        
        # stat.drop_ghost(blocklist)
        # stat.compute_profile( blocklist, bbox, vars_output, 
        #                       outfile=outfile, roughwall=roughwall)
        
        
        stat.extract_profile( xz, bbox, 
                              outfile=outfile, roughwall=roughwall)
        print(f"Profile data saved to {outfile}.\n")

    process_stat(xz_ridge,  bbox_ridge,  f'profile_{loc}_ridge.dat')
    process_stat(xz_valley, bbox_valley, f'profile_{loc}_valley.dat')
    


# =============================================================================
if __name__ == "__main__":

    main()

