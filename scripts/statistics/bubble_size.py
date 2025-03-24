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
import time
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

# from   vista.log         import Logger
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.tools       import get_filelist
#sys.stdout = Logger( os.path.basename(__file__) )

case_folder = '/home/wencan/temp/241030/'

bbox = [-60.0, 100.0, -1.5, 10.0, -11.0, 11.0]

dirs = Directories( case_folder )

parameters = Params( dirs.case_para_file )
roughwall  = parameters.roughwall


with timer('load grid data'):
    grd = GridData( dirs.grid )
    grd.read_grid()
    grd.cell_volume()
    block_list = grd.select_blockgrids( bbox )


if roughwall:
    
    with timer('load wall distance snapshot'):
        wd_snap_file = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]
        wd_snap = Snapshot( wd_snap_file )
        wd_snap.read_snapshot( block_list, var_read=['wd'] )
    
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
    vol1, vol2 = stat.compute_bubble_volume( grd, block_list, cc_df=cc_df,  
                                             roughwall=roughwall,
                                             y_threshold=0.0)

print(f"case {dirs.case_dir} bubble volume (threshold y>0): {vol1:.2f} ({vol2:.2f})")

os.chdir( dirs.pp_bubble )
with open(f"{dirs.pp_bubble}/bubble_size_stat.dat", 'w') as f:
    f.write(f"bubble volume (threshold y>0): {vol1:.2f} ({vol2:.2f})\n")

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
sys.stdout.flush()