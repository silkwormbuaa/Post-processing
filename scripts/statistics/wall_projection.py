#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   get_wd.py
@Time    :   2023/09/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   wall projection script for wavy wall cases
'''

import os
import sys
import pickle
import time

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

import numpy             as     np
import pandas            as     pd

from   vista.statistic   import StatisticData

from   vista.snapshot    import Snapshot

from   vista.directories import Directories

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.plane_analy import save_isolines
from   vista.plane_analy import shift_coordinates
from   vista.plane_analy import periodic_average
from   vista.plane_analy import compute_separation_ratio

from   vista.tools       import get_filelist
from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_wall_projection

from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

bbox = [ -60.0, 108.0, -1.3, 0.01, -11.0, 11.0]

# =============================================================================

dirs = Directories( os.getcwd() )

outpath   = dirs.pp_wall_proj

stat_file = dirs.statistics 
grid_file = dirs.grid
ccfile    = dirs.cc_setup

fric_file = outpath + '/friction_projection.pkl'
pres_file = outpath + '/wall_vars_projection.pkl'

snapshotfile = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]

parametersfile = dirs.case_para_file

# - enter outpath

if not os.path.exists(outpath): 
    os.mkdir( outpath )
    print(f"Created directory {outpath}.\n")

os.chdir(outpath)


if not (os.path.exists(fric_file) and os.path.exists(pres_file)):

    # read in statistics.bin and also grid and cut cell info.
    with timer("read in grid"):
        
        G = GridData( grid_file )
        G.read_grid()
        
        block_list = G.select_blockgrids( bbox, mode='within' )
#        block_list = [885,870,855,840,825,810,795,780]  #upstream blocks

    with timer("read statistics"):
        
        S = StatisticData( stat_file )

        vars = ['u','mu','rho','p','pp']

        with open( stat_file, 'br' ) as f:
            S.read_stat_header( f )
            S.read_stat_body(f, block_list, vars)
        
        S.grid3d = G
        
    # read in wall distance
    with timer("read wall distance field"):
        
        wd_snap = Snapshot( snapshotfile )
        wd_snap.read_snapshot(block_list)

        for num in block_list:
            S.bl[num-1].df['wd'] = wd_snap.snap_data[num-1].df['wd']
            
        print(S.bl[num-1].df)


    # read in cutcell info

    with timer("read cutcell info"):
        
        cc_df = pd.read_csv( ccfile, delimiter=r'\s+')
        cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                    inplace=True)


    with timer("compute friction projection"):
        
        S.friction_projection( block_list, G, cc_df )


    with timer("compute pressure projection"):
        
        S.wall_vars_projection( block_list, G, cc_df )


else:
    
    with timer("read projection data"):
        S = StatisticData( stat_file )
        
        with open( fric_file, 'rb' ) as f:
            S.df_fric = pickle.load( f )
        
        with open( pres_file, 'rb' ) as f:
            S.df_wall = pickle.load( f )

with timer("plotting"):

    parameters = read_case_parameter( parametersfile )
    delta   = float( parameters.get('delta_0') )
    h_ridge = float( parameters.get('H') )
    h_md    = float( parameters.get('H_md') )
    x_imp   = float( parameters.get('x_imp') )
    rho_ref = float( parameters.get('rho_ref') )
    u_ref   = float( parameters.get('u_ref') )
    p_ref   = float( parameters.get('p_ref') )
    period  = int(   parameters.get('period') )
    
    dyn_p   = 0.5*rho_ref*u_ref*u_ref
    
    S.df_fric = shift_coordinates( S.df_fric, delta, h_ridge, h_md, x_imp)
    S.df_wall = shift_coordinates( S.df_wall, delta, h_ridge, h_md, x_imp)
    
    # drop points that before -20.0 delta or after 10.0 delta
    S.df_fric = S.df_fric[ (S.df_fric['xs']>=-20.01) &(S.df_fric['xs']<= 10.01)]
    S.df_wall = S.df_wall[ (S.df_wall['xs']>=-20.01) &(S.df_wall['xs']<= 10.01)]

    xx     = np.array( S.df_fric['xs'] )
    zz     = np.array( S.df_fric['zs'] )
    fric   = np.array( S.df_fric['fric'] )
    p      = np.array( S.df_wall['p'] )
    p_fluc = np.array( S.df_wall['p`'] )

    npx = len( np.unique(xx) )
    npz = len( np.unique(zz) )

    xx     = xx.reshape( npz, npx )
    zz     = zz.reshape( npz, npx )
    fric   = fric.reshape( npz, npx )
    p      = p.reshape( npz, npx )
    p_fluc = p_fluc.reshape( npz, npx )

# --- save original wall projection results

    save_isolines( xx, zz, fric, 1.0, "separationline_xz.pkl")
    
    cbar_levels = np.linspace(-5.0,10.0,31)
    plot_wall_projection( xx, zz, fric/dyn_p*1000.0, 
                          separation="separationline_xz.pkl",
                          filename='fric',
                          cbar_levels=cbar_levels)
    
    cbar_levels = np.linspace( 1.0, 2.25, 26 )
    plot_wall_projection( xx, zz, p/p_ref, 
                          separation="separationline_xz.pkl",
                          filename='pressure',
                          cbar_levels=cbar_levels)
    
    cbar_levels = np.linspace(0.025,0.085,31)
    plot_wall_projection( xx, zz, p_fluc/p_ref,
                          separation="separationline_xz.pkl", 
                          filename='pressure_fluc',
                          cbar_levels=cbar_levels)
    
# --- periodic average results
    
    fric = periodic_average( fric, period )
    
    save_isolines( xx, zz, fric, 1.0, "separationline_xz_periodic.pkl" )

    cbar_levels = np.linspace(-5.0,10.0,31)
    plot_wall_projection( xx, zz, fric/dyn_p*1000.0, 
                          separation="separationline_xz_periodic.pkl",
                          filename='fric_periodic',
                          cbar_levels=cbar_levels)

    cbar_levels = np.linspace( 1.0, 2.25, 26 )
    plot_wall_projection( xx, zz, p/p_ref, 
                          separation="separationline_xz_periodic.pkl",
                          filename='pressure_periodic',
                          cbar_levels=cbar_levels)

# --- output separation area ratio and length ratio distribution

    print(f"separation ratio {compute_separation_ratio(fric):10.5f}.")
#    
#    
    
with timer("save spanwise averaged variable distribution along x"):
    
    fric_mean   = np.mean( fric/dyn_p*1000.0, axis=0 )
    p_mean      = np.mean( p/p_ref, axis=0 )
    p_fluc_mean = np.mean( p_fluc/p_ref, axis=0 )

    df_streamwise = pd.DataFrame(columns=['x','Cf','Cp','p_fluc'])
    df_streamwise['x']  = np.unique( xx )
    df_streamwise['Cf'] = np.array( fric_mean )
    df_streamwise['Cp'] = np.array( p_mean )
    df_streamwise['p_fluc'] = np.array( p_fluc_mean )
    
    with open('streamwise_vars.pkl','wb') as f:
        
        pickle.dump( df_streamwise, f )
    
    df_streamwise.to_string('streamwise_vars.dat',
                            index=False,
                            float_format='%15.7f',
                            justify='left')


# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()      
    
         