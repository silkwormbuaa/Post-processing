#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   compute_correlation_at_shear_layer.py
@Time    :   2025/05/12 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Extract the spanwise u distribution at the pfmax and compute the auto-correlation.

'''

import os
import sys
import time
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.mpi         import MPIenv
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.directories import create_folder

# - build MPI communication environment

mpi = MPIenv()
# =============================================================================

casefolder = '/home/wencan/temp/220927/'

# =============================================================================

dirs      = Directories( casefolder )

params    = Params( dirs.case_para_file )
grd       = None
snapfiles = None
len_prb   = 0
dz        = 0.0
u_stat    = None
p_stat    = None

print(f"Rank {mpi.rank} is working on {dirs.case_dir}.")
sys.stdout.flush()

# - root processor read in grid

if mpi.is_root:
    
    print(f"I am root, now at {dirs.case_dir}.")
    sys.stdout.flush()
    
    snapfiles = get_filelist( dirs.snp_dir, key='snapshot.bin' )
    
    with timer('load grid data'):
        grd = GridData( dirs.grid )
        grd.read_grid()
        
        loc = [params.loc_sl[0]*params.delta_0 + params.x_imp, params.loc_sl[1]*params.delta_0]
        blocklist, indx_probe = grd.select_probed_blockgrids( 'Z', loc )
        for bl_num in blocklist:
            len_prb += grd.g[bl_num-1].nz
        
        dz = grd.g[blocklist[0]-1].dz[0]
            
    print(f"Root reports: snapshots number: {len(snapfiles)}, probe length: {len_prb}, min length: {dz}.")
    
    with timer('load statistic data'):
        stat = StatisticData( dirs.statistics )
        stat.read_statistic( blocklist, vars_in=['u','p'] )
        df_stat = stat.get_probed_df( blocklist, grd, indx_probe, 'Z' )
        u_stat  = df_stat['u'].mean()
        p_stat  = df_stat['p'].mean()
    
    sys.stdout.flush()
    
# - broadcast grid data to all processors

mpi.comm.barrier()
grd       = mpi.comm.bcast( grd,       root=0 )
snapfiles = mpi.comm.bcast( snapfiles, root=0 )
len_prb   = mpi.comm.bcast( len_prb,   root=0 )
dz        = mpi.comm.bcast( dz,        root=0 )
u_stat    = mpi.comm.bcast( u_stat,    root=0 )
p_stat    = mpi.comm.bcast( p_stat,    root=0 )

n_snaps    = len(snapfiles)
u_corrs    = np.zeros( (n_snaps,len_prb), dtype=float )
p_corrs    = np.zeros( (n_snaps,len_prb), dtype=float )
snap_steps = np.zeros( n_snaps,           dtype=int   )
snap_times = np.zeros( n_snaps,           dtype=float )
std_u_t    = np.zeros( n_snaps,           dtype=float ) # std of u at an instantaneous time
std_p_t    = np.zeros( n_snaps,           dtype=float ) # std of p at an instantaneous time
u_avgs     = np.zeros( n_snaps,           dtype=float ) # spanwise average u
p_avgs     = np.zeros( n_snaps,           dtype=float ) # spanwise average p

# - compute the correlation of velocity and pressure in z direction at shear layer

def compute_correlation_at_shear_layer( snapfile ):
    
    index = snapfiles.index( snapfile )
    loc   = [params.loc_sl[0]*params.delta_0 + params.x_imp, params.loc_sl[1]*params.delta_0]
    blocklist, indx_probe = grd.select_probed_blockgrids( 'Z', loc )

    snap = Snapshot( snapfile )
    snap.read_snapshot( blocklist, var_read=['u','p'] )

    df_snap   = snap.get_probed_df( blocklist, grd, indx_probe, 'Z' )

    u_corr = get_auto_correlation( np.array(df_snap['u'])-df_snap['u'].mean() )
    p_corr = get_auto_correlation( np.array(df_snap['p'])-df_snap['p'].mean() )
    
    # plt.figure(figsize=(8,6))
    # plt.plot(df_snap['z'], df_snap['u']-df_snap['u'].mean(), label='u')
    # plt.savefig(f'u_{index}.png', dpi=300)
    # plt.close()

    snap_steps[index]   = snap.itstep
    snap_times[index]   = snap.itime
    u_corrs   [index,:] = u_corr
    p_corrs   [index,:] = p_corr
    
    # plt.figure(figsize=(8,6))
    # plt.plot(np.arange(len(u_corr)), u_corr, label='u')
    # plt.savefig(f'u_cor_{index}.png', dpi=300)
    # plt.close()
    
    # compute the std 
    std_u_t[index] = np.std( df_snap['u'] )
    std_p_t[index] = np.std( df_snap['p'] )
    u_avgs[index]  = df_snap['u'].mean()
    p_avgs[index]  = df_snap['p'].mean()
    

def get_auto_correlation( array:np.ndarray ):
    
    n = len(array)
    mean = np.mean(array)
    array = array - mean
    auto_corr = np.correlate(array, array, mode='full')
    
    corr = auto_corr / auto_corr[n-1]
    
    corr = corr[n-1:]
    return corr 
    

# - distribute works

clock = timer("get pressure at pfmax from snapshots:")

if mpi.size == 1:
    print("No workers available. Master should do all tasks.")
    
    for i, snapfile in enumerate(snapfiles):
        compute_correlation_at_shear_layer( snapfile )
        clock.print_progress( i, len(snapfiles), rank=mpi.rank )

else:
    if mpi.is_root:
        mpi.master_distribute( snapfiles )
    else:
        while True:
            task_index = mpi.worker_receive()
            
            if task_index is None: break
            else:
                compute_correlation_at_shear_layer( snapfiles[task_index] )
                clock.print_progress( task_index, len(snapfiles), rank=mpi.rank )
                
# - gather results

mpi.barrier()

snap_steps = mpi.comm.reduce( snap_steps, root=0, op=mpi.MPI.SUM )
snap_times = mpi.comm.reduce( snap_times, root=0, op=mpi.MPI.SUM )
u_corrs    = mpi.comm.reduce( u_corrs,    root=0, op=mpi.MPI.SUM )
p_corrs    = mpi.comm.reduce( p_corrs,    root=0, op=mpi.MPI.SUM )
std_u_t    = mpi.comm.reduce( std_u_t,    root=0, op=mpi.MPI.SUM )
std_p_t    = mpi.comm.reduce( std_p_t,    root=0, op=mpi.MPI.SUM )
u_avgs     = mpi.comm.reduce( u_avgs,     root=0, op=mpi.MPI.SUM )
p_avgs     = mpi.comm.reduce( p_avgs,     root=0, op=mpi.MPI.SUM )

# - write out the results

if mpi.is_root:
    
    os.chdir( create_folder( dirs.pp_snp_sl_prb ))
    
    u_corrs = np.array(u_corrs).mean(axis=0)
    p_corrs = np.array(p_corrs).mean(axis=0)
    
    lags    = np.arange(0,len(u_corrs))*dz/params.delta_0

    with open('correlations_in_shear_layer.dat','w') as f:
        f.write("wave_length".rjust(15)+"u_corr".rjust(15)+"p_corr".rjust(15)+"\n")
        for i in range(len(u_corrs)):
            f.write(f"{lags[i]:15.6f}{u_corrs[i]:15.6f}{p_corrs[i]:15.6f}\n")
    
    plt.figure(figsize=(8,6))
    plt.plot(lags, u_corrs, label='u_corr')
    plt.plot(lags, p_corrs, label='p_corr')
    plt.xlabel(r'wave length ($\delta_0$)')
    plt.ylabel('correlation')
    plt.title('Auto-correlation in shear layer')
    plt.legend()
    plt.savefig('auto_correlation.png', dpi=300)
    plt.close()
    
    first_u_neg = np.where(u_corrs < 0)[0][0]
    L_u = np.trapz(u_corrs[:first_u_neg], lags[:first_u_neg])
    first_p_neg = np.where(p_corrs < 0)[0][0]
    L_p = np.trapz(p_corrs[:first_p_neg], lags[:first_p_neg])
    
    with open('integral_length.dat','w') as f:
        f.write("L_u  :".rjust(15) + f"{L_u:15.6f}\n")
        f.write("L_p  :".rjust(15) + f"{L_p:15.6f}\n")
        
    # std
    
    # average of std of u and p (spanwise locally)
    std_u = np.mean(std_u_t)
    std_p = np.mean(std_p_t)
    
    # std of spanwise averaged u and p (large scale)
    std_u_avg = np.std(u_avgs)
    std_p_avg = np.std(p_avgs)
    
    u_avg = np.mean(u_avgs)
    p_avg = np.mean(p_avgs)
    
    with open('std.dat','w') as f:
        f.write("std_u  :".rjust(15) + f"{std_u:15.6f}\n")
        f.write("std_p  :".rjust(15) + f"{std_p:15.6f}\n")
        f.write("std_u_avg  :".rjust(15) + f"{std_u_avg:15.6f}\n")
        f.write("std_p_avg  :".rjust(15) + f"{std_p_avg:15.6f}\n")
        f.write("u_avg  :".rjust(15) + f"{u_avg:15.6f}" + f"{u_stat:15.6f}\n")
        f.write("p_avg  :".rjust(15) + f"{p_avg:15.6f}" + f"{p_stat:15.6f}\n")
        
    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}.")
    sys.stdout.flush()
