#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size.py
@Time    :   2024/04/22 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   compute the separation bubble size considering the cutcells
'''


import os
import sys
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.log         import Logger
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.tools       import read_case_parameter
from   vista.timer       import timer
sys.stdout = Logger( os.path.basename(__file__) )

bbox = [-60.0, 90.0, -1.5, 10.0, -11.0, 11.0]

dirs = Directories( os.getcwd() )

parameters = read_case_parameter( dirs.case_para_file )
roughwall  = True if parameters.get('roughwall').lower() == 'true' else False


with timer('load grid data'):
    grd = GridData( dirs.grid )
    grd.read_grid()
    grd.cell_volume()
    block_list = grd.select_blockgrids( bbox)


if roughwall:
    
    with timer('load wall distance snapshot'):
        wd_snap_file = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]
        wd_snap = Snapshot( wd_snap_file )
        wd_snap.read_snapshot( var_read=['wd'] )
    
    with timer('load cutcell info'):
        cc_df = pd.read_csv( dirs.cc_setup, delimiter=r'\s+' )
        cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                    inplace=True)
else: cc_df = None
    
    
with timer('load statistics.bin'):
    
    stat = StatisticData( dirs.statistics )
    
    with open( dirs.statistics, 'rb' ) as f:
        
        stat.read_stat_header( f )
        vars = ['u']
        stat.read_stat_body( f, block_list, vars )
        
        if roughwall:
            stat.assign_wall_dist( block_list, wd_snap )


with timer('compute bubble volume'):
    vol = stat.compute_bubble_volume( grd, block_list, cc_df=cc_df,  roughwall=roughwall )

print(f"case {dirs.case_dir} bubble volume: {vol}")
