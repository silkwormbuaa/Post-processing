#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   read_statbin.py
@Time    :   2023/02/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   read statistic binary data.
           ! Set the headers for cutcell_setup first
'''

import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

import pandas            as     pd
import numpy             as     np

from   vista.statistic   import StatisticData
from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.tools       import get_filelist

# =============================================================================
# option zone
# =============================================================================

bbox = [-57, -49.625, -1.2576, 45.0, -11.0, 11.0]

# =============================================================================

resultspath = os.getcwd()

outpath  = resultspath + '/profile'
datafile = resultspath + '/statistics.bin'
gridfile = resultspath + '/inca_grid.bin'
ccfile   = resultspath + '/cutcells_setup.dat'

snapshotfile = get_filelist(resultspath.split('/results')[0] +'/wall_dist',
                            key='snapshot.bin')[0]

parametersfile = resultspath.split('/results')[0] + '/case_parameters'

outfile  = 'profile_mean.dat'


# - enter outpath

if not os.path.exists(outpath): 
    os.mkdir( outpath )
    print(f"Created directory {outpath}.\n")

os.chdir(outpath)


# - read in grid file

with timer("read in grid"):
    
    G = GridData( gridfile )
    G.verbose = False

    G.read_grid()
    
    # given a rectangular region, get a list of blocks in the region
    block_list = G.select_blockgrids( bbox )


# - read in wall distance data

with timer("read wall distance field"):
    
    wd_snap = Snapshot( snapshotfile )
    wd_snap.read_snapshot( block_list )


# - read in cut cell data and assign vol_fra

with timer("read in cut cell info "):

    cc_df = pd.read_csv( ccfile, delimiter = r'\s+' )

    cc_df.drop( columns=['nx',  'ny',  'nz'
                        ,'fax0','fay0','faz0'
                        ,'fax1','fay1','faz1']
                        , inplace=True )
    
    for num in block_list:

        # dataframe slice for a certain block
        temp_df = cc_df[ cc_df['block_number'] == num ]
        
        wall_dist = np.array( wd_snap.snap_data[num-1][5]['wd'] )
        
        # block number starts from 1, but python list index
        # starts from 0
        G.g[num-1].assign_vol_fra( temp_df, wall_dist )


# - read in statistics data

with timer("read block statistics data "):

    S = StatisticData( datafile )
    
    with open( datafile, 'br' ) as f:   
            
        S.read_stat_header( f )
        vars = ['u','v','w','p','rho','T','uu','vv','ww','uv']
        
        # only blocks in the list will be filled data chunk
        S.read_stat_body( f, block_list, vars )
    
    
with timer("match grid and drop ghost cells"):   
     
    S.match_grid( G, block_list )
    S.drop_ghost( G, block_list )


with timer("compute profile"):
    
    S.compute_profile( block_list, bbox, vars, outfile=outfile )
        