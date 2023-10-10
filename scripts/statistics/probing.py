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

import pandas            as     pd
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.statistic   import StatisticData

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.tools       import read_case_parameter
from   vista.tools       import get_filelist

from   vista.line        import LineData



# =============================================================================

outfolder = '/probing'
probe_type = 'Y'
loc = (-53.6,0.0)

# =============================================================================

datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'

outpath  = datapath + outfolder
parametersfile = datapath.split('/results')[0] + '/case_parameters'

# - read in case parameters

parameters = read_case_parameter( parametersfile )
delta   = float( parameters.get('delta_0') )
casecode =  str( parameters.get('casecode') )

# - read in grid info

G = GridData( gridfile )
G.read_grid()


# - enter outpath

if not os.path.exists(outpath): 
    os.mkdir( outpath )
    print(f"Created directory {outpath}.\n")

os.chdir(outpath)


# - do probing 

block_list, indx_probe = G.select_probed_blockgrids( probe_type, loc )

# - read statistics data file

with timer(" read selected blocks from statistics.bin"):
    
    with open(datafile,'br') as f:
        
        S = StatisticData(datafile)
        S.read_stat_header( f )
        vars = ['u','v','w','uu','vv','ww','uv']
        S.read_stat_body( f, block_list, vars )

        

# - do probing

with timer(" Get probed dataframe"):
    
    df_stat = S.get_probed_df( block_list, G, indx_probe, probe_type )
    
    outfile = 'probing_'+probe_type
    
    df_stat.to_string( outfile, 
                       index=False, 
                       float_format='%15.7f',
                       justify='left' )
        