#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   compute_profile_local.py
@Time    :   2026/01/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Do spanwise periodic average first, then compute the profile at a local point.
'''


import os
import sys
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.params      import Params
from   vista.directories import create_folder


def main():
   
    case_dirs = ['smooth_adiabatic','250821','250710']
    
    loc_norm = -13
    
    for case in case_dirs:
        
        compute_profile_local( '/home/wencan/temp/' + case, loc_norm )
    
    
def compute_profile_local( case_folder, loc_norm ): 
   
    dirs         = Directories( case_folder )
    params       = Params( dirs.case_para_file )

    grid         = GridData( dirs.grid )
    grid.read_grid()
    grid.cell_volume()
    
    vars_read   = ['u','v','w','p','pp','rho','mu','T','uu','vv','ww','uv']
    vars_output = ['u','v','w','p','pp','rho','mu','T','uu','vv','ww','uv']

    outfolder = dirs.pp_statistics + '/profile_' + str(int(loc_norm))
    os.chdir( create_folder( outfolder ) )

    z_valley  = 0.5*params.D
    roughwall = params.roughwall
    
    loc     = params.x_imp + loc_norm*params.delta_0
    xz_DL   = np.array( [loc, 0.001         ] )
    xz_CL   = np.array( [loc, z_valley+0.001] )
    
    bbox_DL = [loc-0.1,loc+0.1, 0.0, 100, -20, 20]
    bbox_CL = [loc-0.1,loc+0.1, 0.0, 100, -20, 20]
    
# =============================================================================

    def process_stat(xz, bbox, outfile):

        blocklist_full = grid.select_blockgrids(bbox, mode='overlap')
        print( blocklist_full )
        
        
        blocklist,_    = grid.select_probed_blockgrids('Y',xz, bbox=bbox, bbox_mode='overlap')
        
        stat_file      = dirs.statistics
        stat           = StatisticData(stat_file)
        stat.grid3d    = grid
        stat.read_statistic(blocklist_full,vars_in=vars_read)
        

        for num in blocklist:
            grid.g[num-1].assign_vol_fra()
        
        stat.match_grid( blocklist_full, grid, add_to_df=True )
                
        # periodic averaging
        stat.spanwise_periodic_average( blocklist_full, vars_output, params.D )
        
        stat.extract_profile( xz, bbox, 
                              outfile=outfile, roughwall=roughwall)
        print(f"Profile data saved to {outfile}.\n")

    process_stat(xz_DL, bbox_DL, f'profile_{str(loc_norm)}_DL.dat')
    process_stat(xz_CL, bbox_CL, f'profile_{str(loc_norm)}_CL.dat')
    


# =============================================================================
if __name__ == "__main__":

    main()

