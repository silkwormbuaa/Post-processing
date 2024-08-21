#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   shock_tracking3d.py
@Time    :   2024/08/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import pandas            as     pd
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.tools       import find_indices
from   vista.tools       import get_filelist

# =============================================================================

y0 = 10.4
xrange = [-10.0, 5.0]

snapshotpath = '/media/wencanwu/Seagate Expansion Drive1/temp/231124/snapshots'
gridfile = '/media/wencanwu/Seagate Expansion Drive1/temp/231124/results/inca_grid.bin'

# =============================================================================

# - get the snapshot file list

snap_files = get_filelist( snapshotpath, 'snapshot.bin' )

# - read grid file

grd = GridData( gridfile )
grd.read_grid()
blocklist = grd.select_blockgrids([xrange[0], xrange[1], y0, y0, -10.4, 10.4])

print(f"These blocks contains the shock-detecting line:{blocklist}.")
print("==================\n")

# - initialize a list to store the x location of the shock

times   = list()
x_shock = list()

# - loop over the snapshot files

for i, snap_file in enumerate(snap_files):
    
    # - read snapshot file

    snap = Snapshot( snap_file )
    snap.verbose = True
    snap.grid3d = grd
    snap.read_snapshot(block_list=blocklist)
    snap.compute_gradients(block_list=blocklist)

    fig, ax = plt.subplots()

    # - initialize a pandas dataframe to store the grad_rho on the probe line

    prbdf = pd.DataFrame( columns=['x','z','grad_rho'] )

    # - loop over the blocks to extract the grad_rho on the probe line
    
    for bn in blocklist:
        
        if bn not in snap.bl_nums:
            raise ValueError(f"block {bn} is not in the snapshot")
        
        else:
            
            bl = snap.snap_data[ snap.bl_nums.index(bn) ]
            g = grd.g[bn-1]        
            
            grad_rho = np.array( bl.df['grad_rho'] ).reshape( bl.npz, bl.npy, bl.npx)

            # find the indices of the probe in y direction
            j,_ = find_indices( g.py, y0 )
            
            plane = np.ravel(grad_rho[3:-3,j,3:-3])
            
            gx,gz = np.meshgrid( g.gx[3:-3], g.gz[3:-3] )
            
            gz = np.ravel(gz)
            gx = np.ravel(gx)
            
            new_data = pd.DataFrame( {'z':gz,'x':gx, 'grad_rho':plane} )

            prbdf = pd.concat( [prbdf, new_data], ignore_index=True )        

    # drop outside data and sort the dataframe
    
    prbdf.drop( prbdf[(prbdf['x'] < xrange[0]) | (prbdf['x'] > xrange[1])].index, 
                inplace=True )
    prbdf = prbdf.sort_values(by=['z','x'])
    prbdf = prbdf.reset_index(drop=True)

    #
    npz = len(np.unique(prbdf['z']))
    npx = len(np.unique(prbdf['x']))
    xx = np.array(prbdf['x']).reshape(npz,npx)
    zz = np.array(prbdf['z']).reshape(npz,npx)
    grad_rho = np.array(prbdf['grad_rho']).reshape(npz,npx)
    
    ax.contourf( xx, zz, grad_rho, cmap='RdBu_r', levels=31)
#        plt.savefig(f'snap_{snap.itstep}.png')

    plt.show()    
    plt.close()

    times.append( snap.itime )        
    
    print(f"{i+1}/{len(snap_files)} is done.")
    print("------------------\n")
