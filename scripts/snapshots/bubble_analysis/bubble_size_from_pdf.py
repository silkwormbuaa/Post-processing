#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size_from_pdf.py
@Time    :   2024/05/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Script of computing the bubble size from the separation pdf.

             Need inca_grid.bin, cutcells_setup.dat, snapshot_container.pkl and wall_dist/ in the same directory.
'''


import os
import sys
import pickle
import numpy             as     np
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import get_filelist

# =============================================================================

case_dir = '/home/wencan/temp/smooth_adiabatic/'
var_out  = ['u','n_sep','pdf_sep','wd']

# =============================================================================

dirs      = Directories( case_dir )
params    = Params( dirs.case_para_file )
roughwall = params.roughwall

# load the grid data

grd       = GridData( dirs.grid )
grd.read_grid()
grd.cell_volume()

# load the snapshot container of separation PDF

os.chdir( dirs.pp_bubble )
with open(f"{dirs.pp_bubble}/snapshot_container.pkl",'rb') as f:
    snapshot_container: Snapshot = pickle.load(f)

# load the wall distance field if roughwall 

if roughwall:
    wd_snap_file = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]
    with timer("load wall distance snapshot data"):
        wd_snap = Snapshot( wd_snap_file )
        wd_snap.read_snapshot( var_read=['wd'] )
    
    snapshot_container.copy_var_from( wd_snap, ['wd'] )

else:
    
    for bl in snapshot_container.snap_data:
        bl_num       = bl.num
        g            = grd.g[bl_num-1]
        X, Y, Z      = np.meshgrid( g.gx, g.gy, g.gz, indexing='ij' )
        bl.df['wd']  = Y.T.flatten()

if roughwall:
    cc_df = pd.read_csv( dirs.cc_setup, delimiter=r'\s+' )
    cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                inplace=True)
else: cc_df = None


# compute the bubble size by the definition of 50% pdf

with timer("compute bubble size"):
    
    sep1 = snapshot_container.compute_bubble_volume_pdf( grd, cc_df, roughwall=roughwall, opt=1 )
    sep2 = snapshot_container.compute_bubble_volume_pdf( grd, cc_df, roughwall=roughwall, opt=2 )

    print(f"bubble size by 50% pdf: {sep1:.2f}")
    print(f"bubble size PDF:        {sep2:.2f}")
    
    with open('bubble_size_pdf.dat','w') as f:
        f.write(f"bubble size by 50% pdf: {sep1:.2f}\n")
        f.write(f"bubble size PDF:        {sep2:.2f}")
    
# write the PDF file into tecplot szplt format.

snapshot_container.grid3d = grd
snapshot_container.write_vtm( 'pdf_sep.vtm', vars=var_out, buff=2 )

#snapshot_container.write_szplt( 'pdf_sep.szplt', buff=2 )
