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

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.statistic   import StatisticData
from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.tools       import get_filelist
from   vista.line        import LineData


# =============================================================================

outfolder = '/spectra'
probe_type = 'Z'
loc = (-42.5,0.52)

# =============================================================================

dirs = Directories( os.getcwd() )

datafile = dirs.statistics 
gridfile = dirs.grid
snapfiles = get_filelist( dirs.snp_dir, 'snapshot.bin' ) 
outpath  = dirs.pos_dir + outfolder

# - read in grid info

G = GridData( gridfile )
G.read_grid()

# - enter outpath

create_folder( outpath )
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

        df_stat = S.get_probed_df( block_list, G, indx_probe, probe_type )


# - read snapshot data file

with timer(" read selected blocks from snapshot.bin"):
    
    for snapfile in snapfiles:
 
        snap3d = Snapshot( snapfile )
        snap3d.read_snapshot( block_list )
        
    # - get probed dataframe
        
        df_snap = snap3d.get_probed_df( block_list, G, indx_probe, probe_type)
        
        df_fluc = pd.DataFrame(columns=['z','u`','v`','w`'])
        
        df_fluc['z']  = df_stat['z']
        df_fluc['u`'] = df_snap['u'] - df_stat['u']
        df_fluc['v`'] = df_snap['v'] - df_stat['v']
        df_fluc['w`'] = df_snap['w'] - df_stat['w']

        
        line = LineData( df=df_fluc )
        
        line.fft_kz( )
        
#        print(line.df)
        
        os.chdir( snapfile.split('/snapshot.bin')[0] )

        line.df.to_string("spectra.dat",
                          index=False,
                          float_format='%15.7f',
                          justify='left')
        
        print(f"Finish {snapfile}")

    
    
    
        