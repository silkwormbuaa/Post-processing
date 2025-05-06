#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   read_statbin.py
@Time    :   2025/04/05
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Integerate the tke in the incipient interaction region
'''

import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

import pandas            as     pd
import numpy             as     np

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.colors      import colors          as    clr
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.tools       import get_filelist

# =============================================================================
# option zone
# =============================================================================

def main():
    
    case_dirs = ['smooth_mid',      '231124','241030','241018',
                 'smooth_adiabatic','221014','220926','220825','220927','221221',
                 '240210',          '240211']
    
    # case_dirs = ['smooth_adiabatic', '220927']

    outfile = '/home/wencan/temp/DataPost/midRe/tke_incip/tke_incip_influence_region.dat'
    with open( outfile, 'w' ) as f:
        f.write( "case     massflow       tke/mf    uu/mf    momen/mf    tke     uu     momen \n")
        
        for case in case_dirs:
            
            mass, tke, uu, momen = compute_tke_incip( '/home/wencan/temp/' + case )
            f.write( f"{case}  {mass:10.5f} {tke:10.5f} {uu:10.5f} {momen:10.5f} {tke*mass:15.5f} {uu*mass:15.5f} {momen*mass:15.5f}\n" )

# =============================================================================

def compute_tke_incip( case_dir ):
    
    print(f"\nStart processing {case_dir}.\n")
        
    dirs = Directories( case_dir )

    outpath  = dirs.pp_profile_incip
    datafile = dirs.statistics
    gridfile = dirs.grid
    ccfile   = dirs.cc_setup

    # - read parameters to check if it is a rough wall case

    params     = Params( dirs.case_para_file )
    roughwall  = params.roughwall
    loc        = (params.x_incip - 2.0)*params.delta_0 + params.x_imp

    if roughwall:
        snapshotfile = get_filelist(dirs.wall_dist,key='snapshot.bin')[0]

    # - enter outpath

    os.chdir( create_folder(outpath) )

# =============================================================================
        
    print(f"Processing extracting profile at location: {loc}.\n")
    
    bbox    = [ loc-2.5, loc+2.5, -1.2576, 8.3,  -11.0, 11.0]  # bounding box

    # - read in grid file

    with timer("read in grid"):
        
        G = GridData( gridfile )
        G.read_grid()
        G.cell_volume()
        
        # given a rectangular region, get a list of blocks in the region
        block_list = G.select_blockgrids( bbox )

    # - read in statistics data

    with timer("read block statistics data "):

        S = StatisticData( datafile )
                
        vars = ['u','v','w','rho','T','uu','vv','ww','uv']
        
        # only blocks in the list will be filled data chunk

        S.read_statistic( block_list, vars )
        S.compute_vars( block_list, ['RS','mach'] )

    # - assign vol_fra to grid, then match G to S

    with timer("Assign vol_fra to G"):
        
        if roughwall:

            # - read in wall distance data
            with timer("read wall distance field"):
                
                wd_snap = Snapshot( snapshotfile )
                wd_snap.read_snapshot( block_list, var_read=['wd'] )

            # - read in cut cell data and assign vol_fra

            with timer("read in cut cell info "):

                cc_df = pd.read_csv( ccfile, delimiter = r'\s+' )

                cc_df.drop( columns=['fax0','faz0'
                                    ,'fax1','faz1', 'processor']
                                    , inplace=True )
                
                for num in block_list:

                    # dataframe slice for a certain block
                    temp_df = cc_df[ cc_df['block_number'] == num ]
                    
                    wall_dist = np.array( wd_snap.snap_data[num-1].df['wd'] )
                    
                    # block number starts from 1, but python list index
                    # starts from 0
                    G.g[num-1].assign_vol_fra( df=temp_df, wall_dist=wall_dist )
                    
        else:
            
            for num in block_list:
                G.g[num-1].assign_vol_fra()

        # match grid and pass vol_fra from G to S
        S.match_grid( block_list, G, add_to_df=True )
        
    # - compute tke

    with timer("compute tke"):
        
        tke      = S.integrate_vol_var( block_list, G, var='tke',  type='massflow', bbox=bbox )
        massflow = S.integrate_vol_var( block_list, G, var=None,   type='massflow', bbox=bbox )
        uu       = S.integrate_vol_var( block_list, G, var='u`u`', type='massflow', bbox=bbox )
        momentum = S.integrate_vol_var( block_list, G, var='u',    type='massflow', bbox=bbox )
        print(clr.fg.green,f"tke = {tke:.5f}.\n"    ,clr.reset)
        print(clr.fg.green,f"massflow = {massflow:.5f}.\n",clr.reset)
        print(clr.fg.green,f"tke/massflow = {tke/massflow:.5f}.\n",clr.reset)

    print(f"Finished computing tke for {params.casecode} at location {loc}.\n")
    
    return massflow, tke/massflow, uu/massflow, momentum/massflow


if __name__ == "__main__":
        
    main()