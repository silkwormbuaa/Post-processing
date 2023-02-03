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

sys.path.append('..')

from   vista.statistic   import StatisticData

from   vista.grid        import GridData


datapath = '/home/wencanwu/my_simulation/temp/results/statistics.bin'
gridpath = '/home/wencanwu/my_simulation/temp/results/inca_grid_3D.bin'
outpath  = '/home/wencanwu/my_simulation/temp/results/'
outfile  = 'header.dat'

#A = StatisticData(datapath)
#
#os.chdir(outpath)
#
#with open(datapath,'br') as f:   
#        
#    A.read_stat_header(f)
#
#    A.read_stat_body(f)

G = GridData( gridpath )

os.chdir( outpath )

with open( gridpath, 'rb' ) as f:
    
    G.verbose = True
    
    G.read_grid_header( f )
    
    G.read_grid_body( f )
