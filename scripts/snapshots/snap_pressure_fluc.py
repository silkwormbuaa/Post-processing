#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   snap_pressure_fluc.py
@Time    :   2024/10/08 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Compute instantaneous pressure fluctuation and visualise it.
'''


import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.statistic   import StatisticData
from   vista.tools       import get_filelist

# =============================================================================
# option
# =============================================================================

bbox      = [-30.0, 999.0, -1.0, 31.0, -999.0, 999.0]
vars_in   = ['p']
vars_out  = ['p','p_fluc']

casedir   = '/home/wencan/temp/smooth_mid'

# =============================================================================

dirs      = Directories( casedir )
snaps_dir = dirs.snp_dir + '/snapshot_01013869'
snapfiles = get_filelist( snaps_dir, 'snapshot.bin' )

grid3d = GridData( dirs.grid )
grid3d.read_grid()

blocklist = grid3d.select_blockgrids( bbox, mode='within' )

stat = StatisticData( dirs.statistics )

stat.read_statistic( blocklist, vars_in )

for i, snapfile in enumerate(snapfiles):
    
    snap = Snapshot( snapfile )
    snap.read_snapshot( var_read=vars_in )

    for bl_num in blocklist:
        
        snapblk = snap.snap_data[snap.bl_nums.index(bl_num)]
        statblk = stat.bl[stat.bl_nums.index(bl_num)]
        snapblk.df['p_fluc'] = snapblk.df['p'] - statblk.df['p']

    os.chdir('/home/wencan/temp/smooth_mid/test/')
    
    snap.write_szplt( 'snapshot_01013869.szplt',vars_out,block_list=blocklist ) 
    print(f'Finished processing {snapfile}')

