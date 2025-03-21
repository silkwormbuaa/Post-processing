#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   wall_projection_sw.py
@Time    :   2023/10/07 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Wall projection script for smooth wall.
'''


import os
import sys
import pickle
import time
import numpy             as     np
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.plane_analy import save_isolines
from   vista.plane_analy import shift_coordinates
from   vista.plane_analy import compute_separation_ratio
from   vista.plot_style  import plot_wall_projection
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

bbox = [ -120.0, 108.0, 0, 0.001, -11.0, 11.0]

# =============================================================================

dirs = Directories( os.getcwd() )

outpath   = dirs.pp_wall_proj
stat_file = dirs.statistics
grid_file = dirs.grid

pres_file = outpath + '/wall_vars_projection.pkl'

# - enter outpath

create_folder( outpath )
os.chdir(outpath)

# - check if the wall projection file exists

if not os.path.exists(pres_file) :

    # read in statistics.bin and also grid and cut cell info.
    with timer("read in grid"):
        
        G = GridData( grid_file )
        G.read_grid()
        
        block_list = G.select_blockgrids( bbox, mode='overlap' )
#        block_list = [885,870,855,840,825,810,795,780]  #upstream blocks

    with timer("read statistics"):
        
        S = StatisticData( stat_file )

        vars = ['u','mu','rho','p','pp']

        with open( stat_file, 'br' ) as f:
            S.read_stat_header( f )
            S.read_stat_body(f, block_list, vars)
        
        S.grid3d = G
        
    with timer("extract wall variables"):
        
        S.extract_wall_vars_sw( block_list, G )      

else:
    
    with timer("read projection data"):
        
        S = StatisticData( stat_file )
        
        with open( pres_file, 'rb' ) as f:
            S.df_wall = pickle.load( f )

with timer("plotting"):

    params   = Params( dirs.case_para_file )
    delta    = params.delta_0
    h_ridge  = params.H
    h_md     = params.H_md
    x_imp    = params.x_imp
    rho_ref  = params.rho_ref
    u_ref    = params.u_ref
    p_ref    = params.p_ref
    period   = params.n_period
    casecode = params.casecode
    
    dyn_p   = 0.5*rho_ref*u_ref*u_ref
    
    S.df_wall = shift_coordinates( S.df_wall, delta, h_ridge, h_md, x_imp)
    
    # drop points that before -20.0 delta or after 10.0 delta
    S.df_wall = S.df_wall[ (S.df_wall['xs']>=-20.01) &(S.df_wall['xs']<= 10.01)]

    xx     = np.array( S.df_wall['xs'] )
    zz     = np.array( S.df_wall['zs'] )
    fric   = np.array( S.df_wall['fric'] )
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
    cbar_ticks  = np.linspace(-5.0,10.0,4)
    plot_wall_projection( xx, zz, fric/dyn_p*1000.0, 
                          separation="separationline_xz.pkl",
                          filename='fric',
                          cbar_label=r'$C_f\times 10^3$',
                          label=casecode,
                          cbar_ticks=cbar_ticks,
                          cbar_levels=cbar_levels)
    
    cbar_levels = np.linspace( 1.0, 2.5, 31 )
    cbar_ticks  = np.linspace( 1.0, 2.5, 7 )
    plot_wall_projection( xx, zz, p/p_ref, 
                          separation="separationline_xz.pkl",
                          cbar_label=r'$p/p_{\infty}$',
                          label=casecode,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          filename='pressure' )
    
    cbar_levels = np.linspace(0.020,0.090,36)
    cbar_ticks  = np.linspace(0.020,0.090,8)
    plot_wall_projection( xx, zz, p_fluc/p_ref, 
                          separation="separationline_xz.pkl",
                          cbar_label=r"$\sqrt{\langle p'p'\rangle}/p_{\infty}$",
                          label=casecode,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          filename='pressure_fluc' )


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
    