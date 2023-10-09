#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   spectra.py
@Time    :   2023/10/09 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   From probed fluctuation to get spectra
'''


import os

import sys

import pandas            as     pd
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.statistic   import StatisticData

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.tools       import read_case_parameter
from   vista.tools       import get_filelist

from   vista.line        import LineData


# =============================================================================

outfolder = '/spectra'
probe_type = 'Z'
loc = (-42.25,0.52)

# =============================================================================

datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'

snappath = datapath.split('/results')[0]+'/snapshots'
snapfile = get_filelist( snappath, 'snapshot.bin' )[0]

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
        S.compute_vars( block_list, ['RS'] )


# - read snapshot data file

with timer(" read selected blocks from snapshot.bin"):
    
    with open(snapfile,'rb') as f:
        
        snap3d = Snapshot( snapfile )
        snap3d.read_snapshot( block_list )
        

# - get probed dataframe

with timer("Get probed dataframe "):
    
    pd.set_option('display.max_rows', None)
    
    df_stat = S.get_probed_df( block_list, G, indx_probe, probe_type )
    
    print(df_stat)
    
    df_snap = snap3d.get_probed_df( block_list, G, indx_probe, probe_type)
    
    print(df_snap)
    
    df_fluc = pd.DataFrame(columns=['u`','v`','w`'])
    
    df_fluc['z']  = df_stat['z']
    df_fluc['u`'] = df_snap['u'] - df_stat['u']
    df_fluc['v`'] = df_snap['v'] - df_stat['v']
    df_fluc['w`'] = df_snap['w'] - df_stat['w']

    print(df_fluc)     
    
    line = LineData( df=df_fluc )
    
    line.fft('z','u`')
    
    print(line.df)
    
    fig,ax = plt.subplots()
    ax.loglog(line.df['k_z'],line.df['E_u`'])
    plt.show()
    
    
    
    
        