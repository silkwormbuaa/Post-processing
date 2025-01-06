#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   wall_projection.py
@Time    :   2024/09/23 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   project instantaneous wall friction for rough wall cases.
'''


import gc
import os
import sys
import pickle
import time
import numpy             as     np
import pandas            as     pd
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.material    import get_visc
from   vista.plot_style  import plot_wall_projection
from   vista.directories import create_folder
from   vista.plane_analy import save_isolines
from   vista.plane_analy import shift_coordinates


# - build MPI communication environment

comm    = MPI.COMM_WORLD
rank    = comm.Get_rank()
n_procs = comm.Get_size()

# =============================================================================

casepath = '/home/wencan/temp/241030' 
bbox     = [ -30.0, 108.0, -1.3, 0.0, -11.0, 11.0]     # [xmin, xmax, ymin, ymax, zmin, zmax]

# =============================================================================
# preparation
# =============================================================================

dirs       = Directories( casepath )
grid_file  = dirs.grid
ccfile     = dirs.cc_setup
outpath_cf = dirs.pp_snp_fricprj + '/figs_cf'
outpath_p  = dirs.pp_snp_fricprj + '/figs_p'
outpath_pf = dirs.pp_snp_fricprj + '/figs_pf'

# --- broadcast the parameters

params     = None
roughwall  = True
snapfiles  = None
block_list = None
cc_df      = None
df_stat    = None
grid3d     = GridData()
wd_snap    = Snapshot()

if rank == 0:
    
    create_folder( outpath_p )
    create_folder( outpath_cf )
    create_folder( outpath_pf )
    
    params    = Params( dirs.case_para_file )
    roughwall = params.roughwall
    
    snapfiles = get_filelist( dirs.snp_dir, 'snapshot.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot files.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    block_list = grid3d.select_blockgrids( bbox, mode='within' )
    
    if roughwall:
        
        wd_file = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]
        wd_snap = Snapshot( wd_file )
        wd_snap.read_snapshot( block_list, var_read=['wd'] )

        cc_df   = pd.read_csv( ccfile, delimiter=r'\s+')
        cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], inplace=True)

    # read in the time-averaged wall projection dataframe
    
    with open( dirs.pp_wall_proj+'/wall_vars_projection.pkl','rb') as f:
        df_stat = pickle.load(f)


params     = comm.bcast( params,    root=0 )
roughwall  = comm.bcast( roughwall, root=0 )
snapfiles  = comm.bcast( snapfiles, root=0 )
block_list = comm.bcast( block_list,root=0 )
cc_df      = comm.bcast( cc_df,     root=0 )
df_stat    = comm.bcast( df_stat,   root=0 )
grid3d     = comm.bcast( grid3d,    root=0 )
wd_snap    = comm.bcast( wd_snap,   root=0 )


Re_ref     = params.Re_ref
visc_law   = params.visc_law

# --- distribute tasks

n_snaps    = len( snapfiles )
i_s, i_e   = distribute_mpi_work(n_snaps, n_procs, rank)
snapfiles  = snapfiles[i_s:i_e]

print(f"I am processor {rank:05d}, I take {len(snapfiles):5d} tasks.")
sys.stdout.flush()

comm.barrier()

# --- loop over the snapshots

clock = timer("show cf")

for i, snap_file in enumerate(snapfiles):
    
    snap3d = Snapshot( snap_file )
    snap3d.grid3d = grid3d
    snap3d.read_snapshot( block_list=block_list, var_read=['u','p','rho','T'] )

    itstep  = snap3d.itstep
    itime   = snap3d.itime
    
    if roughwall:
        snap3d.copy_var_from( wd_snap, ['wd'], block_list )

    for bl in snap3d.snap_data:

        if bl.num in block_list:
            bl.df['mu'] = get_visc( np.array(bl.df['T']), Re_ref, law=visc_law )

# ----- wall projection

    if roughwall:
        df_wall  = snap3d.friction_projection( block_list, grid3d, cc_df )
        df_wall2 = snap3d.wall_vars_projection( block_list, grid3d, cc_df )
        df_wall['p'] = df_wall2['p']
    else:
        df_wall = snap3d.extract_wall_vars_sw( block_list, grid3d )   

# ----- visualization

    delta    = params.delta_0
    h_ridge  = params.H
    h_md     = params.H_md
    x_imp    = params.x_imp
    rho_ref  = params.rho_ref
    u_ref    = params.u_ref
    p_ref    = params.p_ref
    
    df_stat  = shift_coordinates( df_stat, delta, h_ridge, h_md, x_imp)
    df_wall  = shift_coordinates( df_wall, delta, h_ridge, h_md, x_imp)
    dyn_p    = 0.5*rho_ref*u_ref*u_ref
    
    # drop points that before -14.5 delta or after 10.0 delta
    df_stat  = df_stat[ (df_stat['xs']>=-14.5) &(df_stat['xs']<= 10.01)]    
    df_wall  = df_wall[ (df_wall['xs']>=-14.5) &(df_wall['xs']<= 10.01)]

    xx       = np.array( df_wall['xs'] )
    zz       = np.array( df_wall['zs'] )
    fric     = np.array( df_wall['fric'] )
    p        = np.array( df_wall['p'] )
    p_mean   = np.array( df_stat['p'] )
    p_fluc   = p - p_mean

    npx      = len( np.unique(xx) )
    npz      = len( np.unique(zz) )
    

    xx       = xx.reshape( npz, npx )
    zz       = zz.reshape( npz, npx )
    fric     = fric.reshape( npz, npx )
    p        = p.reshape( npz, npx )
    p_fluc   = p_fluc.reshape( npz, npx )
    
# --- save original wall projection results

    os.chdir( outpath_cf )

    save_isolines( xx, zz, fric, 0.0, f"sep_line_{itstep:08d}.pkl")

    figname = f'cf_{itstep:08d}.png'
    cbar_levels = np.linspace(-6.0,6.0,41)
    cbar_ticks  = np.linspace(-6.0,6.0,5)
    plot_wall_projection( xx, zz, fric/dyn_p*1000.0, 
                          separation=f"sep_line_{itstep:08d}.pkl",
                          filename=figname,
                          col_map='RdBu_r',
                          cbar_label=r'$C_f\times 10^3$',
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          use_norm=True,
                          title=f"t = {itime:8.2f} ms")
    os.remove(f"sep_line_{itstep:08d}.pkl")

    os.chdir( outpath_p )
    
    cbar_levels = np.linspace(0.8,2.4,41)
    cbar_ticks  = np.linspace(0.8,2.4,5)
    plot_wall_projection( xx, zz, p/p_ref,
                          filename = f'p_{itstep:08d}.png',
                          col_map='RdBu_r',
                          extend='both',
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          title=f"t = {itime:8.2f} ms")

    os.chdir( outpath_pf )

    cbar_levels = np.linspace(-0.4,0.4,41)
    cbar_ticks  = np.linspace(-0.4,0.4,5)
    plot_wall_projection( xx, zz, p_fluc/p_ref,
                          filename = f'p_fluc_{itstep:08d}.png',
                          col_map='RdBu_r',
                          extend='both',
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          use_norm=True,
                          title=f"t = {itime:8.2f} ms")
    

# - print the progress

    del snap3d, df_wall
    gc.collect()
    
    progress = (i+1)/len(snapfiles)
    print(f"Rank:{rank:05d},{i+1}/{len(snapfiles)} is done. " + clock.remainder(progress))
    print("------------------\n")
    sys.stdout.flush()


# print out the time finishing the job

comm.barrier()

if rank == 0:

    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush()  