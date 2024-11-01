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
from   vista.math_opr    import find_parabola_max

# =============================================================================

y0 = 31.2
xrange = [17,32]    

# y = 10.4, xrange = [-8,10] for 231124
#           xrange = [8,24]  for smooth_mid
# y = 31.2, xrange = [17,32] for 231124
#           xrange = [30,45] for smooth_mid


inputpath    = '/home/wencan/temp/231124/'

outputpath   = inputpath + 'postprocess/snapshots/shock_tracking/2d/'
snapshotpath = inputpath + 'snapshots/'
gridfile     = inputpath + 'results/inca_grid.bin'

tolerance  = 3.0     # max distance between max_grad_rho and mean shock location

half_width = 2.0    # half width of the subdomain to search for the shock front

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

times         = list()
x_shocks      = list()
grad_rho_maxs = list()
x_last_shock  = None

# - loop over the snapshot files

os.chdir( outputpath )
clock = timer("Tracking the shock location")

for i, snap_file in enumerate(snap_files):
    
    # - read snapshot file

    snap = Snapshot( snap_file )
    snap.grid3d = grd
    snap.read_snapshot()
    snap.compute_gradients()

    fig, ax = plt.subplots(figsize=(8,6))

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

            if prbdf.empty: prbdf = new_data
            else:           prbdf = pd.concat( [prbdf, new_data], ignore_index=True )        

    # - drop outside data and sort the dataframe
    
    prbdf.drop( prbdf[(prbdf['x'] < xrange[0]) | (prbdf['x'] > xrange[1])].index, 
                inplace=True )
    prbdf = prbdf.sort_values(by='x')
    prbdf = prbdf.reset_index(drop=True)

    # - find the maximum of the gradient
    
    idmax = prbdf['grad_rho'].idxmax()
    
    p1 = [prbdf['x'][idmax-1], prbdf['grad_rho'][idmax-1]]
    p2 = [prbdf['x'][idmax],   prbdf['grad_rho'][idmax]]
    p3 = [prbdf['x'][idmax+1], prbdf['grad_rho'][idmax+1]]
    
    x_shock, grad_rho_max = find_parabola_max(p1,p2,p3)
    
    # x_shock = prbdf.iloc[idmax]['x']
    # grad_rho_max = prbdf.iloc[idmax]['grad_rho']

# -----------------------------------------------------------------------------
# - check if the maximum is at the boundary
    
    if idmax == 0 or idmax == len(prbdf)-1:
        print("The maximum is at the boundary! Please have a check!")

# - check if x_shock is 'tolerance' away from x_last_shock

    if x_last_shock is None: x_last_shock = x_shock

    if abs(x_shock - x_last_shock) > tolerance:
        print("Warning: the shock front is not continuous! Special treatment will be applied!")

        indx_s,_ = find_indices( np.array(prbdf['x']), x_last_shock-half_width )
        indx_e,_ = find_indices( np.array(prbdf['x']), x_last_shock+half_width )
        sub_df   = prbdf.iloc[indx_s:indx_e].copy()

        idmax = sub_df['grad_rho'].idxmax()

        p1 = [sub_df['x'][idmax-1], sub_df['grad_rho'][idmax-1]]
        p2 = [sub_df['x'][idmax],   sub_df['grad_rho'][idmax]]
        p3 = [sub_df['x'][idmax+1], sub_df['grad_rho'][idmax+1]]
        
        x_shock, grad_rho_max = find_parabola_max(p1,p2,p3)

        # x_shock      = sub_df['x'][idmax]
        # grad_rho_max = sub_df['grad_rho'][idmax]
        
# -----------------------------------------------------------------------------
    
    # - plot the grad_rho on the probe line
    
    ax.plot( prbdf['x'], prbdf['grad_rho'], label='grad_rho' )
    ax.plot( x_shock, grad_rho_max, 'ro' )
    
    ax.set_xlim( xrange )
    ax.set_ylim( [-0.05, 0.8] )
    ax.set_title(f"t= {snap.itime:8.2f} ms")
    
    plt.savefig(f'snap_Z_{snap.itstep:08d}.png')
    
    plt.close()
    
    # - output

    x_last_shock = x_shock
    times.append( snap.itime )        
    x_shocks.append( x_shock )
    grad_rho_maxs.append( grad_rho_max )
    
    progress = (i+1)/len(snap_files)
    print(f"{i+1}/{len(snap_files)} is done. " + clock.remainder(progress))
    print("------------------\n")

result = pd.DataFrame( {'time':times, 'x_shock':x_shocks, 'grad_rho_max':grad_rho_maxs} )
result.to_csv('shock_tracking2d.csv', index=False, sep=' ')