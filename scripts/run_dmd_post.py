# -*- coding: utf-8 -*-
'''
@File    :   run_dmd_post.py
@Time    :   2023/06/05 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os

import sys

import numpy             as     np

source_dir = os.path.dirname( os.path.dirname( os.path.realpath(__file__) ))
sys.path.append( source_dir )

from   vista.timer       import timer

from   vista.snapshot    import Snapshot

from   vista.dmdmodes    import DMDModes

from   vista.tools       import get_filelist

# read in one snapshot file to get grid vectors

snap_dir = os.getcwd()

with timer('\nGet snapshots file and grid vector'):
    
    snap_files = get_filelist( snap_dir + '/snapshots' )
    
    testfile = snap_files[0]
    
    snapshot_temp = Snapshot( testfile )
    
    x,y = snapshot_temp.get_grid_vectors(buff=3)

    print(f"shape of x:{np.shape(x)}, from [{x[:3]}...{x[-3:]}]")
    print(f"shape of y:{np.shape(y)}, from [{y[:3]}...{y[-3:]}]")


# match data and grids:


