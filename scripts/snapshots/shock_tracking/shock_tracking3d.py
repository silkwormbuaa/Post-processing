#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   shock_tracking3d.py
@Time    :   2024/08/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Slice at y= 2*delta0 and track the shock front. 
             All use special treatment.(only search in the subdomain around last shock)
'''


import os
import sys
import pickle
import pandas            as     pd
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.mpi         import MPIenv
from   vista.timer       import timer
from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.params      import Params
from   vista.directories import Directories
from   vista.tools       import find_indices
from   vista.tools       import distribute_mpi_work
from   vista.directories import create_folder
from   vista.math_opr    import find_parabola_max


def main():
    
    mpi = MPIenv()
    
    case_dir = '/home/wencan/temp/smooth_mid/'
    
    dirs     = Directories( case_dir )
    params   = Params( dirs.case_para_file )

    grd       = None
    snapfiles = None 

    if mpi.is_root:
        
        create_folder( dirs.pp_shock )
        create_folder( dirs.pp_shock + '/group1' )
        create_folder( dirs.pp_shock + '/group2' )
        grd = GridData( dirs.grid )
        grd.read_grid()
        snapfiles = dirs.snap3d_files

    grd       = mpi.comm.bcast(grd, root=0)
    snapfiles = mpi.comm.bcast(snapfiles, root=0)
    params = Params( dirs.case_para_file )

    os.chdir( dirs.pp_shock + '/group1' )
    #mpi_shock_tracking( mpi, grd, snapfiles, params.shock_range_1, 'shock_tracking1.pkl' )
    static_alloc_shock_tracking( mpi, grd, snapfiles, params.shock_range_1, 'shock_tracking1.pkl' )
    
    mpi.barrier()
    
    os.chdir( dirs.pp_shock + '/group2' )
    #mpi_shock_tracking( mpi, grd, snapfiles, params.shock_range_2, 'shock_tracking2.pkl' )    
    static_alloc_shock_tracking( mpi, grd, snapfiles, params.shock_range_2, 'shock_tracking2.pkl' )

# ----------------------------------------------------------------------
# >>> static_alloc_shock_tracking                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/03/21  - created
#
# Desc
#
# ----------------------------------------------------------------------

def static_alloc_shock_tracking( mpi, grd, snapfiles, ranges, outfile ):
    
    if mpi.is_root: 
        create_folder( './figs/'); sys.stdout.flush()

    mpi.barrier()
    
    clock = timer("Tracking the shock front")
    
    blocklist = grd.select_blockgrids([ranges[1], ranges[2], ranges[0], ranges[0], -10.4, 10.4])
    
    # - initialize a list to store the x location of the shock
    times                 = list()
    shocklines            = list()
    shockline             = None
    
    # - distribute the tasks among all processors
    
    i_s, i_e = distribute_mpi_work( len(snapfiles), mpi.size, mpi.rank )
    snapfiles = snapfiles[i_s:i_e]
    
    for i, snapfile in enumerate(snapfiles):

        itime, shockline = snap_shock_tracking( 
                           snapfile, 
                           grd, 
                           blocklist, 
                           ranges, 
                           shockline )
        times.append( itime )
        shocklines.append( shockline )
        
        clock.print_progress( i, len(snapfiles), mpi.rank )
    
    gather_time = mpi.comm.gather( times, root=0 )
    gather_line = mpi.comm.gather( shocklines, root=0 )

    if mpi.is_root:
        
        times = [ t for sublist in gather_time for t in sublist ]
        shocklines = [ l for sublist in gather_line for l in sublist ]
        
        # sort based on time
        times, shocklines = zip(*sorted(zip(times, shocklines)))
        
        with open(outfile, 'wb') as f:
            pickle.dump( times, f )
            pickle.dump( shocklines, f )
        
        print(f"Task {i} completed.")
        sys.stdout.flush()
        

# =============================================================================
# - do shock tracking in paralel

# !!! realize dynamic task allocation is not applicable here, since the shock
# !!! tracking is dependent on the last shock location.

def mpi_shock_tracking( mpi, grd, snapfiles, ranges, outfile ):
    
    blocklist = grd.select_blockgrids([ranges[1], ranges[2], ranges[0], ranges[0], -10.4, 10.4])

    # - initialize a list to store the x location of the shock
    times                 = list()
    shocklines            = list()
    shockline             = None
    
    if mpi.is_root: 
        create_folder( './figs/'); sys.stdout.flush()

    mpi.barrier()
    
    clock = timer("Tracking the shock front")
    
    if mpi.size == 1:
    
        print("No worker available. Master should do all tasks.")
        sys.stdout.flush()
        
        for i, snapfile in enumerate(snapfiles):
            
            itime, shockline = snap_shock_tracking( 
                               snapfile, 
                               grd, 
                               blocklist, 
                               ranges, 
                               shockline )
            times.append( itime )
            shocklines.append( shockline )

            clock.print_progress( i, len(snapfiles) )

        gather_time = mpi.comm.gather( times, root=0 )
        gather_line = mpi.comm.gather( shocklines, root=0 )

    else:
        
        if mpi.is_root:
            mpi.master_distribute( snapfiles )
            
        else:
            while True:
                task_id = mpi.worker_receive()
                if task_id is None: break
                else:
                    itime, shockline = snap_shock_tracking( 
                                       snapfiles[task_id], 
                                       grd, 
                                       blocklist, 
                                       ranges, 
                                       shockline )
                    times.append( itime )
                    shocklines.append( shockline )

                    clock.print_progress( task_id, len(snapfiles) )

        mpi.barrier()
        
        gather_time = mpi.comm.gather( times, root=0 )
        gather_line = mpi.comm.gather( shocklines, root=0 )

    if mpi.is_root:
        
        times = [ t for sublist in gather_time for t in sublist ]
        shocklines = [ l for sublist in gather_line for l in sublist ]
        
        # sort based on time
        times, shocklines = zip(*sorted(zip(times, shocklines)))
        
        with open(outfile, 'wb') as f:
            pickle.dump( times, f )
            pickle.dump( shocklines, f )


# =============================================================================
# - do shock tracking on a single snapshot

def snap_shock_tracking( snap_file, grid, blocklist, ranges, x_last_shock,
                         tolerance=3.0, half_width=2.0 ):
    
    # - read snapshot file

    snap = Snapshot( snap_file )
    snap.verbose = False
    snap.grid3d  = grid
    snap.read_snapshot(block_list=blocklist)
    snap.compute_gradients(block_list=blocklist)

# - initialize a pandas dataframe list to store all blocks' df on the probe line

    prbdf_list = list()

    # - loop over the blocks to extract the grad_rho on the probe line
    
    for bn in blocklist:
        
        if bn not in snap.bl_nums:
            raise ValueError(f"block {bn} is not in the snapshot")
        
        else:
            
            bl = snap.snap_data[ snap.bl_nums.index(bn) ]
            g = grid.g[bn-1]        
            
            grad_rho = np.array( bl.df['grad_rho'] ).reshape( bl.npz, bl.npy, bl.npx)

            # find the indices of the probe in y direction
            j,_ = find_indices( g.py, ranges[0] )
            
            plane = np.ravel(grad_rho[3:-3,j,3:-3])
            
            gx,gz = np.meshgrid( g.gx[3:-3], g.gz[3:-3] )
            
            gz = np.ravel(gz)
            gx = np.ravel(gx)
            
            new_data = pd.DataFrame( {'z':gz,'x':gx, 'grad_rho':plane} )
            prbdf_list.append( new_data )
    
    # - concatenate the dataframes
    
    prbdf = pd.concat( [df for df in prbdf_list], ignore_index=True )        

    # - drop outside data and sort the dataframe
    
    prbdf.drop( prbdf[(prbdf['x'] < ranges[1]) | (prbdf['x'] > ranges[2])].index, 
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
    
    idmax  = grad_rho.argmax(axis=1)
    
    x_temp = xx[0,idmax]

# - check if any element of x_shock is 'tolerance' away from x_last_shock
# ------------------------------------------------------------------------------    
    
    if x_last_shock is None: 
        
        x_last_shock = x_temp
        if np.any( np.abs(x_temp - np.mean(x_last_shock))  > tolerance ):
            print("Warning: the shock front is not continuous! Special treatment will be applied.\n")
            x_last_shock = np.array([np.mean(x_temp)]*len(x_last_shock))

    else: x_last_shock = x_last_shock['x'].values
    
    indx_s    = np.array([find_indices(xx[0,:], x-half_width)[0] for x in x_last_shock])
    indx_e    = indx_s + int( 2*half_width//abs(xx[0,1]-xx[0,0]) )
    
    indices      = np.array([np.arange(s,e) for s,e in zip(indx_s,indx_e)])
    sub_grad_rho = np.take_along_axis(grad_rho, indices, axis=1)
    sub_xx       = np.take_along_axis(xx, indices, axis=1)
    
    idmax        = sub_grad_rho.argmax(axis=1) # array stores the index of max grad_rho
    x_shock      = np.zeros(npz)
    
    for j in range(len(idmax)):

        p1 = [ sub_xx[j,idmax[j]-1], sub_grad_rho[j,idmax[j]-1] ]
        p2 = [ sub_xx[j,idmax[j]  ], sub_grad_rho[j,idmax[j]  ] ]
        p3 = [ sub_xx[j,idmax[j]+1], sub_grad_rho[j,idmax[j]+1] ]
        
        x_shock[j], _ = find_parabola_max(p1,p2,p3)        
    
    z_shock   = zz[:,0]
    shockline = pd.DataFrame( {'x':x_shock, 'z':z_shock} )
    
    # - plot the schlieren plane and the shock line
    
    fig, ax = plt.subplots(figsize=(10,6))
    
    clevels = np.linspace(0,0.5,51)
    
    ax.contourf( xx, zz, grad_rho, cmap='Greys', levels=clevels, extend='both' )
    ax.plot( sub_xx[:,-1], z_shock, 'b', ls=':' )
    ax.plot( sub_xx[:,0 ], z_shock, 'b', ls=':' )
    ax.plot( x_shock, z_shock, 'r', ls=':' )
    ax.set_title(f"t ={snap.itime:8.2f} ms")

    plt.savefig(f'./figs/snap_{snap.itstep:08d}.png')
    plt.close()
    
    return snap.itime, shockline
# =============================================================================

if __name__ == '__main__':
    
    main()
