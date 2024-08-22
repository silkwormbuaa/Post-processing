#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   shock_tracking.py
@Time    :   2024/08/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   tracking the separation shock location in the Z snapshots
'''

import os
import sys
import pandas            as     pd
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.timer       import timer
from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.tools       import find_indices
from   vista.tools       import get_filelist

# =============================================================================

y0 = 10.4
xrange = [-10.0, 5.0]

snapshotpath = 'path/to/snapshots'
gridfile = 'path/to/inca_grid.bin'
outpath = 'path/to/output'

# =============================================================================

# - get the snapshot file list

snap_files = get_filelist( snapshotpath, 'snapshot_Z_001.bin' )

# - read grid file

grd = GridData( gridfile )
grd.read_grid()
blocklist = grd.select_blockgrids([xrange[0], xrange[1], y0, y0, 0.01, 0.10])

print(f"These blocks contains the shock-detecting line:{blocklist}.")
print("==================\n")

# - initialize a list to store the x location of the shock

times   = list()
x_shock = list()

# - loop over the snapshot files

os.chdir( outpath )
clock = timer("Tracking the shock location")

for i, snap_file in enumerate(snap_files):
    
    # - read snapshot file

    snap = Snapshot( snap_file )
    snap.grid3d = grd
    snap.read_snapshot()
    snap.compute_gradients()

    fig, ax = plt.subplots()

    # - initialize a pandas dataframe to store the grad_rho on the probe line

    prbdf = pd.DataFrame( columns=['x','grad_rho'] )

    # - loop over the blocks to extract the grad_rho on the probe line
    
    for bn in blocklist:
        
        if bn not in snap.bl_nums:
            raise ValueError(f"block {bn} is not in the snapshot")
        
        else:
            
            bl = snap.snap_data[ snap.bl_nums.index(bn) ]
            g = grd.g[bn-1]        
            
            grad_rho = np.array( bl.df['grad_rho'] ).reshape( bl.npy, bl.npx)

            # find the indices of the probe in y direction
            j,_ = find_indices( g.py, y0 )
            
            line = grad_rho[j,:]
            
            new_data = pd.DataFrame( {'x':g.gx[3:-3], 'grad_rho':line[3:-3]} )

            prbdf = pd.concat( [prbdf, new_data], ignore_index=True )        

    # drop outside data and sort the dataframe
    
    prbdf.drop( prbdf[(prbdf['x'] < xrange[0]) | (prbdf['x'] > xrange[1])].index, 
                inplace=True )
    prbdf = prbdf.sort_values(by='x')
    prbdf = prbdf.reset_index(drop=True)

    # find the maximum of the gradient
    
    idmax = prbdf['grad_rho'].idxmax()

    if idmax == 0 or idmax == len(prbdf)-1:
        
        print("The maximum is at the boundary! Please have a check!")
        ax.plot( prbdf['x'], prbdf['grad_rho'], label='grad_rho' )
        plt.savefig(f'snap_{snap.itstep}.png')
    
    plt.close()

    times.append( snap.itime )        
    x_shock.append( prbdf.iloc[idmax]['x'] )
    
    progress = (i+1)/len(snap_files)
    print(f"{i+1}/{len(snap_files)} is done. " + clock.remainder(progress))
    print("------------------\n")

result = pd.DataFrame( {'time':times, 'x_shock':x_shock} )
result.to_csv('shock_tracking.csv', index=False)