#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   get_block_grid_dic.py
@Time    :   2023/08/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Use this scripts to check the order of blocks
             (statistics.bin and grid.bin, theirs should equal)
             Found the order of statistics blocks increases
             one by one. Therefore, now it is unnecessary to 
             make a dictionary between blocks and grids. 
'''


import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from vista.statistic     import StatisticData

from vista.grid          import GridData

working_dir = '/home/wencanwu/my_simulation/temp/220927_lowRe/results'

os.chdir( working_dir )

S = StatisticData( 'statistics.bin' )

S.verbose = True

with open( 'statistics.bin', 'rb' ) as f:
    
    S.read_stat_header( f )
    
    S.read_stat_body( f, [] )
    
num_bl = len( S.bl )

print(f"number of blocks is {num_bl}.\n")

list_bl_nums = []

for bl in S.bl:
    
    list_bl_nums.append( bl.num )

print( list_bl_nums )
    
