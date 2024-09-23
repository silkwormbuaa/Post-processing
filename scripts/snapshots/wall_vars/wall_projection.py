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
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.tools       import read_case_parameter
from   vista.material    import get_visc
from   vista.plot_style  import plot_wall_projection
from   vista.directories import create_folder
from   vista.plane_analy import save_isolines
from   vista.plane_analy import shift_coordinates
from   vista.plane_analy import compute_separation_ratio


# - build MPI communication environment

comm    = MPI.COMM_WORLD
rank    = comm.Get_rank()
n_procs = comm.Get_size()

# =============================================================================

bbox     = [ -30.0, 108.0, -1.3, 0.01, -11.0, 11.0]     # [xmin, xmax, ymin, ymax, zmin, zmax]
casepath = '/home/wencanwu/test/220927' 

# =============================================================================
# preparation
# =============================================================================

dirs      = Directories( casepath )
grid_file = dirs.grid
ccfile    = dirs.cc_setup
outpath   = dirs.pp_snp_fricprj
wd_file   = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]

snapfiles  = None
block_list = None
cc_df      = None
params     = None
p_dyn      = None
grid3d     = GridData()
wd_snap    = Snapshot()

if rank == 0:
    
    create_folder( outpath )
    
    params  = read_case_parameter( dirs.case_para_file )
    u_ref   = float(params.get('u_ref'))
    rho_ref = float(params.get('rho_ref'))
    p_dyn   = 0.5 * rho_ref * u_ref**2
    
    snapfiles = get_filelist( dirs.snp_dir, 'snapshot.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot files.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()
    block_list = grid3d.select_blockgrids( bbox, mode='within' )
    
    wd_snap = Snapshot( wd_file )
    wd_snap.read_snapshot( block_list, var_read=['wd'] )

    cc_df = pd.read_csv( ccfile, delimiter=r'\s+')
    cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], inplace=True)

    
p_dyn      = comm.bcast( p_dyn,      root=0 )
snapfiles  = comm.bcast( snapfiles,  root=0 )
block_list = comm.bcast( block_list, root=0 )
params     = comm.bcast( params,     root=0 )
grid3d     = comm.bcast( grid3d,     root=0 )
wd_snap    = comm.bcast( wd_snap,    root=0 )
cc_df      = comm.bcast( cc_df,      root=0 )

n_snaps   = len( snapfiles )
i_s, i_e  = distribute_mpi_work(n_snaps, n_procs, rank)
snapfiles = snapfiles[i_s:i_e]

print(f"I am processor {rank:05d}, I take {len(snapfiles):5d} tasks.")
sys.stdout.flush()

os.chdir( outpath )
clock = timer("show cf")

# - loop over the snapshots

for i, snap_file in enumerate(snapfiles):
    
    snap3d = Snapshot( snap_file )
    snap3d.grid3d = grid3d
    snap3d.verbose = True
    snap3d.read_snapshot( block_list=block_list, var_read=['u','T'] )

    itstep  = snap3d.itstep
    itime   = snap3d.itime
    df_file = f'df_fric_{itstep:08d}.pkl'
    figname = f'cf_{itstep:08d}.png'
    
    snap3d.copy_var_from( wd_snap, ['wd'] )

    for bl in snap3d.snap_data:

        if bl.num in block_list:
            bl.df['mu'] = get_visc( np.array(bl.df['T']) )

# --- wall projection

    snap3d.friction_projection( block_list, grid3d, cc_df )     

# --- visualization

    delta    = float( params.get('delta_0') )
    h_ridge  = float( params.get('H') )
    h_md     = float( params.get('H_md') )
    x_imp    = float( params.get('x_imp') )
    rho_ref  = float( params.get('rho_ref') )
    u_ref    = float( params.get('u_ref') )
    p_ref    = float( params.get('p_ref') )
    period   = int(   params.get('period') )
    casecode = params.get('tag')
    
    dyn_p   = 0.5*rho_ref*u_ref*u_ref
    
    df_fric = shift_coordinates( snap3d.df_fric, delta, h_ridge, h_md, x_imp)
    
    # drop points that before -20.0 delta or after 10.0 delta
    df_fric = df_fric[ (df_fric['xs']>=-20.01) &(df_fric['xs']<= 10.01)]

    xx     = np.array( df_fric['xs'] )
    zz     = np.array( df_fric['zs'] )
    fric   = np.array( df_fric['fric'] )

    npx    = len( np.unique(xx) )
    npz    = len( np.unique(zz) )

    xx     = xx.reshape( npz, npx )
    zz     = zz.reshape( npz, npx )
    fric   = fric.reshape( npz, npx )

# --- save original wall projection results

    save_isolines( xx, zz, fric, 1.0, "separationline_xz.pkl")
    
    cbar_levels = np.linspace(-5.0,10.0,31)
    cbar_ticks  = np.linspace(-5.0,10.0,4)
    plot_wall_projection( xx, zz, fric/dyn_p*1000.0, 
                          separation="separationline_xz.pkl",
                          filename='fric',
                          cbar_label=r'$C_f\times 10^3$',
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          label=casecode)
    
# --- output separation area ratio and length ratio distribution

    print(f"separation ratio {compute_separation_ratio(fric):10.5f}.")

    fric_mean   = np.mean( fric/dyn_p*1000.0, axis=0 )

    df_streamwise = pd.DataFrame(columns=['x','Cf'])
    df_streamwise['x']  = np.unique( xx )
    df_streamwise['Cf'] = np.array( fric_mean )
    
    df_streamwise.to_string('streamwise_vars.dat',
                            index=False,
                            float_format='%15.7f',
                            justify='left')

# - print the progress
    
    progress = (i+1)/len(snapfiles)
    print(f"Rank:{rank:05d},{i+1}/{len(snapfiles)} is done. " + clock.remainder(progress))
    print("------------------\n")


# print out the time finishing the job

comm.barrier()

if rank == 0:

    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush()  