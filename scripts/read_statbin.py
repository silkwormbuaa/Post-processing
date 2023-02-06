#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   read_statbin.py
@Time    :   2023/02/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   read statistic binary data.
'''

import sys

import os

import re

import pandas            as     pd

import numpy             as     np

sys.path.append('..')

from   vista.statistic   import StatisticData

from   vista.grid        import GridData

from   utils.timer       import timer


datapath = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/statistics.bin'
gridpath = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/inca_grid.bin'
outpath  = '/home/wencanwu/my_simulation/temp/results/'
outfile  = 'header.dat'
ccfile   = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/cutcells_setup.dat'

#%%
'''
A = StatisticData(datapath)

#os.chdir(outpath)

with open(datapath,'br') as f:   
        
    A.read_stat_header(f)

    A.read_stat_body(f,[3,4])
'''
#%% read grid

G = GridData( gridpath )

#os.chdir( outpath )

with open( gridpath, 'rb' ) as f:
    
    G.verbose = True
    
    G.read_grid_header( f )
    
    G.read_grid_body( f )
    
    G.get_sorted_groups()
    
    rect1 = [-57, -1.2576, -49.625, 0]
    
    G.select_blockgrids( rect1 )
    
    block_list = np.array( G.blockgrids_sel ).ravel()

    with timer("read in cut cell info "):

        cc_df = pd.read_csv( ccfile, delimiter = r'\s+' )

        cc_df.drop( columns=['nx',  'ny',  'nz'
                            ,'fax0','fay0','faz0'
                            ,'fax1','fay1','faz1']
                            , inplace=True )
    
    for num in block_list:

        # dataframe slice for a certain block
        temp = cc_df[cc_df['block_number'] == num ]
        
        # block number starts from 1, but python list index
        # starts from 0
        
        G.g[num-1].assign_vol_fra( temp )

