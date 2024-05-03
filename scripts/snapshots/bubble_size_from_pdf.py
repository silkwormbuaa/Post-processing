#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size_from_pdf.py
@Time    :   2024/05/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Script of computing the bubble size from the separation pdf.
'''


import os
import sys
import pickle
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.tools       import get_filelist

# =============================================================================
roughwall = False
# =============================================================================

workpath = os.getcwd()

# load the snapshot container of separation PDF

with open('snapshot_container.pkl','rb') as f:
    snapshot_container: Snapshot = pickle.load(f)


# load the wall distance field if roughwall 

if roughwall:
    wd_snap_file = get_filelist( workpath + '/wall_dist', key='snapshot.bin')[0]
    with timer("load wall distance snapshot data"):
        wd_snap = Snapshot( wd_snap_file )
        wd_snap.read_snapshot( var_read=['wd'] )
        

# load the grid data

grid_file = workpath + '/inca_grid.bin'
cc_file = workpath + '/cutcells_setup.dat'

grd = GridData( grid_file)
grd.read_grid()
grd.cell_volume()

if roughwall:
    cc_df = pd.read_csv( cc_file, delimiter=r'\s+' )
    cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                inplace=True)
else: cc_df = None


# compute the bubble size by the definition of 50% pdf



# write the PDF file into tecplot szplt format.

snapshot_container.drop_ghost()
snapshot_container.write_snapshot_szplt( 'pdf_sep_test.szplt' )



