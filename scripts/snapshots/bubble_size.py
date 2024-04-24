#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size.py
@Time    :   2024/04/22 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   compute the instanteneous separation bubble size considering the cutcells
'''


import os
import sys
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.log         import Logger
from   vista.tools       import get_filelist
from   vista.timer       import timer
sys.stdout = Logger( os.path.basename(__file__) )


dirs = Directories( '/media/wencanwu/Seagate Expansion Drive1/temp/220927' )

grid_file = dirs.grid
ccfile    = dirs.cc_setup
wd_snap_file = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]

snapshotfiles = get_filelist( dirs.snp_dir, key='snapshot.bin')

with timer('load grid data'):
    grd = GridData( grid_file )
    grd.read_grid()
    grd.cell_volume()

with timer('load wall distance snapshot data'):
    
    wd_snap = Snapshot( wd_snap_file )
#    wd_snap.verbose = True
    wd_snap.read_snapshot( var_read=['wd'] )
    
    
with timer("read cutcell info"):
    
    cc_df = pd.read_csv( ccfile, delimiter=r'\s+')
    cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                inplace=True)

for snapshotfile in snapshotfiles:
    
    with timer(f'load snapshot data ...{snapshotfile[-15:]}'):
        
        snap = Snapshot( snapshotfile )
#        snap.verbose = True
        snap.read_snapshot( var_read=['u'] )
        snap.assign_wall_dist( wd_snap)
    
    with timer('compute bubble size'):
        
        vol = snap.compute_bubble_volume( grd, cc_df, roughwall=True )
        
        print(f"snapshot ...{snapshotfile[-25:]} bubble size:{vol:15.4f}")