#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   probing.py
@Time    :   2023/10/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Script of probing data along a line
'''

import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder

# =============================================================================

case_dir   = '/home/wencan/temp/250218'
probe_type = 'X'
loc        = ( 1.31, 0.01 )
outfile    = f'probing_2.dat'

# =============================================================================

dirs = Directories( case_dir )

datafile = dirs.statistics

outpath  = dirs.pp_statistics + '/probing_x'

# - read in case parameters

parameters = Params( dirs.case_para_file )
delta      = parameters.delta_0
casecode   = parameters.casecode

# - read in grid info

G = GridData( dirs.grid )
G.read_grid()

# - enter outpath

os.chdir( create_folder( outpath ) )

# - do probing 

block_list, indx_probe = G.select_probed_blockgrids( probe_type, loc )

# - read statistics data file

with timer(" read selected blocks from statistics.bin"):
    
    stat = StatisticData(dirs.statistics)
    stat.read_statistic( block_list, stat.full_vars )
    
# - do probing

with timer(" Get probed dataframe"):
    
    df_stat = stat.get_probed_df( block_list, G, indx_probe, probe_type )
    
    df_stat.to_string( outfile, 
                       index=False, 
                       float_format='%15.7f',
                       justify='left' )
        