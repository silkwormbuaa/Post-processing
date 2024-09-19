#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   shock_tracking3d.py
@Time    :   2024/08/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Slice at y= 2*delta0 and track the shock front. 
'''


import os
import sys
import pickle
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
xrange = [-8.0, 6.0]

inputpath    = '/media/wencan/Expansion/temp/231124/'
outputpath   = '/home/wencanwu/temp/'

snapshotpath = inputpath + 'snapshots/'
gridfile     = inputpath + 'results/inca_grid.bin'

tolerance = 3.0     # max distance between max_grad_rho and mean shock location

half_width = 1.5    # half width of the subdomain to search for the shock front

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

times      = list()
shocklines = list()
x_last_shock = None

# - loop over the snapshot files

os.chdir( outputpath )
clock = timer("Tracking the shock front")

for i, snap_file in enumerate(snap_files):
    
    # - read snapshot file

    snap = Snapshot( snap_file )
    snap.verbose = False
    snap.grid3d = grd
    snap.read_snapshot(block_list=blocklist)
    snap.compute_gradients(block_list=blocklist)

    fig, ax = plt.subplots(figsize=(10,6))

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

    # - drop outside data and sort the dataframe
    
    prbdf.drop( prbdf[(prbdf['x'] < xrange[0]) | (prbdf['x'] > xrange[1])].index, 
                inplace=True )
    prbdf = prbdf.sort_values(by=['z','x'])
    prbdf = prbdf.reset_index(drop=True)

    # - compose the schlieren plane
    
    npz = len(np.unique(prbdf['z']))
    npx = len(np.unique(prbdf['x']))
    xx = np.array(prbdf['x']).reshape(npz,npx)
    zz = np.array(prbdf['z']).reshape(npz,npx)
    grad_rho = np.array(prbdf['grad_rho']).reshape(npz,npx)
    
    # - tracking the shock front line where max grad_rho is located
    
    idmax = grad_rho.argmax(axis=1)
    
    x_shock   = xx[np.arange(npz),idmax]
    
# - check if any element of x_shock is 'tolerance' away from x_last_shock
# ------------------------------------------------------------------------------    
    
    if x_last_shock is None: x_last_shock = np.mean(x_shock)
    
    if any( abs(x_shock - x_last_shock)  > tolerance ):
        print("Warning: the shock front is not continuous! Special treatment will be applied.\n")
    
        indx_s,_     = find_indices(xx[0,:], x_last_shock-half_width)
        indx_e,_     = find_indices(xx[0,:], x_last_shock+half_width)
        sub_grad_rho = grad_rho[:,indx_s:indx_e]
        sub_xx       = xx[:,indx_s:indx_e]
        sub_zz       = zz[:,indx_s:indx_e]
        
        idmax        = sub_grad_rho.argmax(axis=1)
    
        x_shock   = sub_xx[np.arange(npz),idmax]
        
# ------------------------------------------------------------------------------
    
    z_shock   = zz[:,0]
    shockline = pd.DataFrame( {'x':x_shock, 'z':z_shock} )
    
    x_last_shock = np.mean(x_shock)
    
    # - plot the schlieren plane
    
    clevels = np.linspace(0,0.36,37)
    
    ax.contourf( xx, zz, grad_rho, cmap='Greys', levels=clevels, extend='both' )
    ax.plot( x_shock, z_shock, 'r', ls=':' )
    ax.set_title(f"t ={snap.itime:8.2f} s")

    plt.savefig(f'snap_{snap.itstep:08d}.png')
    plt.close()

    # - output 

    times.append( snap.itime )
    shocklines.append( shockline )
    
    progress = (i+1)/len(snap_files)
    print(f"{i+1}/{len(snap_files)} is done. " + clock.remainder(progress))
    print("------------------\n")

with open('shock_tracking3d.pkl', 'wb') as f:
    pickle.dump( times, f )
    pickle.dump( shocklines, f )
