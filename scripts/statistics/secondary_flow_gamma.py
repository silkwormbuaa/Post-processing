#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   secondary_flow_gamma.py
@Time    :   2024/06/23 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Compute the secondary flow intensity Gamma (Guo et al. 2022)
'''

import os
import sys
import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import griddata

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.tools       import get_filelist

from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

locs_delta = np.linspace(-20,-20,1)
y_lim      = 5.2
outfolder  = '/yz_planes'
compressible = False          # if True, considering density change

# =============================================================================

dirs = Directories( os.getcwd() )

datafile = dirs.statistics
gridfile = dirs.grid
ccfile   = dirs.cc_setup
snapshotfile = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]
outpath  = dirs.pp_statistics + outfolder

# - read in case paramters

params   = Params( dirs.case_para_file )
delta    = params.delta_0
h_ridge  = params.H
h_md     = params.H_md
x_imp    = params.x_imp
p_ref    = params.p_ref
u_ref    = params.u_ref
casecode = params.casecode
n_period = params.n_period
tag      = params.tag

locs = locs_delta*delta + x_imp

# - read in grid info

G = GridData( gridfile )
G.read_grid()

# - read in cutcell info

cc_df = pd.read_csv( ccfile, delimiter=r'\s+')
cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], inplace=True )

# - enter outpath

create_folder( outpath )
os.chdir(outpath)

# - read in statistics data and do slicing

for i,loc in enumerate(locs):
    
    loc_delta = locs_delta[i]
    
    # determine the blocks at the probing location and assign 
    block_list, indx_slic = G.select_sliced_blockgrids( 'X', loc )
    print(f"selected {len(block_list)} blocks at x = {loc}.\n")
    
    # - read in wall distance file and assign volume fraction

    wd_snap = Snapshot( snapshotfile )
    wd_snap.read_snapshot( block_list )
    
    for num in block_list:
    
        # dataframe slice for a certain block
        temp_df = cc_df[ cc_df['block_number'] == num ]
        
        wall_dist = np.array( wd_snap.snap_data[num-1].df['wd'] )
        
        # block number starts from 1, but python list index
        # starts from 0
        G.g[num-1].assign_vol_fra( df=temp_df, wall_dist=wall_dist )

    # - read in statistics data
    
    S = StatisticData( datafile )
    
    with timer("read selected blocks and match grid"):
        
        with open( datafile, 'br' ) as f:
            
            S.read_stat_header( f )
            
            vars = ['v','w','rho']
            
            S.read_stat_body( f, block_list, vars )
            
            S.match_grid( block_list, G )
            
            dfs = S.get_slice_df( block_list, G, indx_slic, 'X' )
            
            dfs.drop( dfs[ dfs['y'] > y_lim ].index, inplace=True )
    
    # - compute secondary flow intensity
    
    with timer("compute secondary flow intensity"):
        
        gamma = 0.0
        
        hy = np.array( dfs['hy'] )
        hz = np.array( dfs['hz'] )
        vol_fra = np.array( dfs['vol_fra'] )
        v   = np.array( dfs['v'] )
        w   = np.array( dfs['w'] )
        rho = np.array( dfs['rho'] )
        
        if not compressible: rho = 1.0
        
        t1 = np.sum( hy*hz*vol_fra * np.sqrt(v*v+w*w) * rho / 507)
        t2 = np.sum( hy*hz*vol_fra * rho)
        
        gamma = t1/t2

    print("=====================================================")
    print(f"Secondary flow intensity at x = {loc} is {gamma:.4f}.\n")
