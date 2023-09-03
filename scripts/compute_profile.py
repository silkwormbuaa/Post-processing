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

from   vista.grid        import GridData

from   vista.timer       import timer

#datapath = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results/'

datapath = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/'

datafile = datapath + 'statistics.bin'
gridfile = datapath + 'inca_grid.bin'

outpath  = datapath

ccfile   = datapath + 'cutcells_setup.dat'

outfile  = 'mean_profile_test.dat'

# - select which wavy wall case
#
#   1 : 1014 case, D/delta = 2
#   2 : 0926 case, D/delta = 1
#   3 : 0825 case, D/delta = 0.5
#   4 : 0927 case, D/delta = 0.25
#   5 : 1221 case, D/delta = 0.125

geo_case = 4

G = GridData( gridfile )

#os.chdir( outpath )

G.verbose = False

## read in whole grid info

# 1. read_grid() : read_grid_header() + read_grid_body()
#   * grid headers: containing what will be read
#   * grid body: every block's grids information
# 2. sorted_group: sort and group block index basd on x,y,z

G.read_grid()

G.get_sorted_groups()


# given a rectangular region, get a list of blocks within the region

bbox = [-57, -49.625, -1.2576, 45.0, -11.0, 11.0]

block_list = G.select_blockgrids( bbox )

with timer("read in cut cell info "):

    cc_df = pd.read_csv( ccfile, delimiter = r'\s+' )

    cc_df.drop( columns=['nx',  'ny',  'nz'
                        ,'fax0','fay0','faz0'
                        ,'fax1','fay1','faz1']
                        , inplace=True )
    
with timer("assign volume fractions "):
    
    for num in block_list:

        # dataframe slice for a certain block
        temp_df = cc_df[cc_df['block_number'] == num ]
        
        # block number starts from 1, but python list index
        # starts from 0
        
        G.g[num-1].assign_vol_fra( temp_df, geo_case )

S = StatisticData( datafile )

with timer("read block statistics data "):
    
    with open( datafile, 'br' ) as f:   
            
        S.read_stat_header( f )
        
        vars = ['u','v','w','p','rho','T','uu','vv','ww','uv']
        
        # only blocks in the list will be filled data chunk
        S.read_stat_body( f, block_list, vars )
    
    
with timer("match grid and drop ghost cells"):   
     
    S.match_grid( G, block_list )
    
    S.drop_ghost( G, block_list )


with timer("compute profile"):
    
    os.chdir(datapath)

    S.compute_profile( block_list, bbox, vars )
        