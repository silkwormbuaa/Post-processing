#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   slice_3dsnapshot.py
@Time    :   2023/08/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os

import sys

import numpy             as     np

import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

snapdir = '/home/wencanwu/my_simulation/temp/220926_lowRe/snapshots/snapshot_00452401/snapshot_block'

os.chdir( snapdir )

""" snapshot3d = Snapshot( 'snapshot.bin' )

snapshot3d.verbose = True

with timer("\nread 3d snapshot and get snapshot struct"):
    
    snapshot3d.get_snapshot_struct()

# write snapshot

with timer("\nwrite 3d snapshot"):
    
    snapshot3d.write_snapshot("snapshot3d.bin")

# get slice of snapshot3d

with timer("\ndo slice"):
        
    snapshot2d = snapshot3d.get_slice( 'Y', 0.0 )
    
    snapshot2d.verbose = True
    
# write snapshots
    
with timer("\nwrite 2d snapshot"):
    
    snapshot2d.write_snapshot("snapshot2d.bin") """

with timer("\nread in new snapshot and show"):
    
    snapshot2d_new = Snapshot("snapshot_W_003.bin")
    
    snapshot2d_new.verbose = True
    
    snapshot2d_new.read_snapshot()
    
    snapshot2d_new.verbose = True
    
    snapshot2d_new.drop_ghost( buff=3 )
    
    snapshot2d_new.assemble_block()
    
    print(snapshot2d_new.df)
    
    x = np.array( snapshot2d_new.df['x'] )
    
    z = np.array( snapshot2d_new.df['z'] )
    
    p = np.array( snapshot2d_new.df['p'] )

    N_x = len(np.unique(x))
    N_z = len(np.unique(z))
    
    x = x.reshape(N_z,N_x)
    z = z.reshape(N_z,N_x)
    p = p.reshape(N_z,N_x)
    
    fig, ax = plt.subplots()
    contour = ax.pcolor(x,z,p)
    ax.set_title('pressure')
    plt.show()
    


