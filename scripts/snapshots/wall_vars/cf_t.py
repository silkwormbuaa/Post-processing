#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   cf_t.py
@Time    :   2025/10/16 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   like wall_projection.py but not plotting, just saving Cf and finally give a Cf_t plot.
'''


import gc
import os
import sys
import pickle
import time
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.mpi         import MPIenv
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.material    import get_visc
from   vista.directories import create_folder
from   vista.plane_analy import save_isolines
from   vista.plane_analy import shift_coordinates

def main():

    mpi = MPIenv()
    
    case_dir = '/home/wencan/temp/smooth_adiabatic'
    
    # be careful with the range of bounding box
    bbox     = [ -30.0, 120.0, -1.3, 0.5, -11.0, 11.0]

    # =============================================================================
    # preparation
    # =============================================================================

    dirs       = Directories( case_dir )
    ccfile     = dirs.cc_setup
    outpath    = dirs.pp_snp_fricprj
    
    if mpi.rank == 0:
        os.chdir( create_folder(outpath) )
        if os.path.exists( 'cf_t.pkl' ):
            read_plot()
            return

    # --- broadcast the parameters

    params     = None
    roughwall  = True
    snapfiles  = None
    block_list = None
    cc_df      = None
    df_stat    = None
    grid3d     = GridData()
    wd_snap    = Snapshot()
    data       = []
    
    if mpi.rank == 0:
        
        create_folder( outpath )
        
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

    params     = mpi.comm.bcast( params,    root=0 )
    roughwall  = mpi.comm.bcast( roughwall, root=0 )
    snapfiles  = mpi.comm.bcast( snapfiles, root=0 )
    block_list = mpi.comm.bcast( block_list,root=0 )
    cc_df      = mpi.comm.bcast( cc_df,     root=0 )
    df_stat    = mpi.comm.bcast( df_stat,   root=0 )
    grid3d     = mpi.comm.bcast( grid3d,    root=0 )
    wd_snap    = mpi.comm.bcast( wd_snap,   root=0 )

    Re_ref     = params.Re_ref
    visc_law   = params.visc_law

    # --- distribute tasks

    n_snaps    = len( snapfiles )
    i_s, i_e   = distribute_mpi_work(n_snaps, mpi.size, mpi.rank)
    snapfiles  = snapfiles[i_s:i_e]

    print(f"I am processor {mpi.rank:05d}, I take {len(snapfiles):5d} tasks.")
    sys.stdout.flush()
    
    mpi.barrier()
    
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
        
        # drop points that before -14.5 delta or after 10.0 delta
        df_stat  = df_stat[ (df_stat['xs']>=-14.5) &(df_stat['xs']<= 10.01)]    
        df_wall  = df_wall[ (df_wall['xs']>=-14.5) &(df_wall['xs']<= 10.01)]
        
        dyn_p    = 0.5*rho_ref*u_ref*u_ref
        
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

        data.append((itstep, itime, 
                     np.mean(fric,axis=0)/dyn_p*1000,
                     np.mean(p,axis=0)/p_ref,
                     np.mean(p_fluc,axis=0)/p_ref ))

    # - print the progress

        del snap3d, df_wall
        gc.collect()
        
        progress = (i+1)/len(snapfiles)
        print(f"Rank:{mpi.rank:05d},{i+1}/{len(snapfiles)} is done. " + clock.remainder(progress))
        print("------------------\n")
        sys.stdout.flush()

    # gather data from all processors
    mpi.barrier()
    all_data = mpi.comm.gather( data, root=0 )
    
    if mpi.rank == 0:
        
        # merge all data
        merged = [item for sub in all_data for item in sub]
        merged.sort( key=lambda x: x[0] )  # sort by itstep
        
        os.chdir( outpath )
        
        itimes = np.array( [ item[1] for item in merged ] )
        cf_t   = np.array( [ item[2] for item in merged ] )
        p_t    = np.array( [ item[3] for item in merged ] )
        pf_t   = np.array( [ item[4] for item in merged ] )
        
        with open('cf_t.pkl','wb') as f:
            pickle.dump( np.unique(xx),f )
            pickle.dump( merged, f )
            
        xx,zz = np.meshgrid( np.unique(xx), itimes )
        
        plot_breathing( xx, zz, cf_t, 
                        cbar_label=r'$C_f\times 10^3$',
                        cbar_levels=np.linspace(-6.0,6.0,81),
                        figname='cf_t',
                        extend='both' )
        
        plot_breathing( xx, zz, p_t, 
                        cbar_label=r'$p/p_{\infty}$',
                        cbar_levels=np.linspace(0.8,2.4,81),
                        figname='p_t',
                        extend='both' )
        
        plot_breathing( xx, zz, pf_t, 
                        cbar_label=r'$p^{\prime}/p_{\infty}$',
                        cbar_levels=np.linspace(-0.4,0.4,81),
                        figname='pf_t',
                        extend='both' )
        

        print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
        sys.stdout.flush()

def read_plot():
    
    print("cf_t.pkl already exists. Just plot it directly.")
    
    with open('cf_t.pkl','rb') as f:
        x_coords   = pickle.load(f)
        merged     = pickle.load(f)
        
    itimes = np.array( [ item[1] for item in merged ] )
    cf_t   = np.array( [ item[2] for item in merged ] )
    p_t    = np.array( [ item[3] for item in merged ] )
    pf_t   = np.array( [ item[4] for item in merged ] )
    
    xx,zz = np.meshgrid( x_coords, itimes )
    
    plot_breathing( xx, zz, cf_t, 
                    cbar_label=r'$C_f\times 10^3$',
                    cbar_levels=np.linspace(-3.5,3.5,71),
                    figname='cf_t',
                    extend='both',
                    u0=True)
    
    plot_breathing( xx, zz, p_t, 
                    cbar_label=r'$p/p_{\infty}$',
                    cbar_levels=np.linspace(0.8,2.4,81),
                    figname='p_t',
                    extend='both' )
    
    plot_breathing( xx, zz, pf_t, 
                    cbar_label=r'$p^{\prime}/p_{\infty}$',
                    cbar_levels=np.linspace(-0.4,0.4,81),
                    figname='pf_t',
                    extend='both' )

def plot_breathing( xx, zz, v, 
                    cbar_label, cbar_levels, 
                    figname,
                    extend='both',
                    u0=False):
    
    fig, ax = plt.subplots( figsize=(15, 50), constrained_layout=True )
    
    cs = ax.contourf( xx, zz, v, levels=cbar_levels, extend=extend, cmap='RdBu_r' )
    if u0:
        u0line = ax.contour( xx, zz, v, levels=[0.0], colors='black', linewidths=0.2 )
    cbar = fig.colorbar( cs, ax=ax, orientation='vertical', pad=0.02 )
    cbar.set_label( cbar_label, fontsize=16 )
    cbar.ax.tick_params( labelsize=14 )
    
    ax.set_xlabel( r"$x/\delta_0$", fontsize=16 )
    ax.set_ylabel( r"$t$ (ms)", fontsize=16 )
    ax.tick_params( axis='both', which='major', labelsize=14 )
    
    plt.savefig( figname + '.png' )
    plt.close()
    

# =============================================================================
if __name__ == "__main__":

    main()

