#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   get_wd.py
@Time    :   2023/09/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import pickle
import numpy             as     np
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.statistic   import StatisticData

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.plane_analy import save_isolines
from   vista.plane_analy import shift_coordinates

from   vista.tools       import get_filelist
from   vista.tools       import read_case_parameter

from   vista.plot_style  import plot_wall_projection

from   vista.lib.form    import phy
from   vista.lib.form    import mth

bbox = [ -60.0, 100.0, -1.3, 0.01, -11.0, 11.0]

resultspath = os.getcwd()

stat_file = resultspath + '/statistics.bin'
grid_file = resultspath + '/inca_grid.bin'
ccfile    = resultspath + '/cutcells_setup.dat'

snapshotfile = get_filelist(resultspath.split('/results')[0] +'/wall_dist',
                            key='snapshot.bin')[0]

parametersfile = resultspath.split('/results')[0] + '/case_parameters'

fric_file = resultspath + '/friction_projection.pkl'
pres_file = resultspath + '/pressure_projection.pkl'

if not (os.path.exists(fric_file) and os.path.exists(pres_file)):

    # read in statistics.bin and also grid and cut cell info.
    with timer("read in grid"):
        G = GridData( grid_file )
        G.read_grid()
        
        block_list = G.select_blockgrids( bbox, mode='within' )


    with timer("read statistics"):
        S = StatisticData( stat_file )

        vars = ['u','mu','p','pp']

        with open( stat_file, 'br' ) as f:
            S.read_stat_header( f )
            S.read_stat_body(f, block_list, vars)
        
        S.grid3d = G
        
        
    # read in wall distance

    wd_snap = Snapshot( snapshotfile )
    wd_snap.read_snapshot(block_list)

    for num in block_list:
        S.bl[num-1].df['wd'] = wd_snap.snap_data[num-1][5]['wd']
        
    print(S.bl[num-1].df)


    # read in cutcell info

    with timer("read cutcell info"):
        
        cc_df = pd.read_csv( ccfile, delimiter=r'\s+')
        cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                    inplace=True)


    with timer("compute friction projection"):
        
        S.friction_projection( block_list, G, cc_df )


    with timer("compute pressure projection"):
        
        S.pressure_projection( block_list, G, cc_df )


else:
    
    with timer("read projection data"):
        S = StatisticData( stat_file )
        
        with open( fric_file, 'rb' ) as f:
            S.df_fric = pickle.load( f )
        
        with open( pres_file, 'rb' ) as f:
            S.df_pres = pickle.load( f )

with timer("plotting"):

    parameters = read_case_parameter( parametersfile )
    delta   = float( parameters.get('delta_0') )
    h_ridge = float( parameters.get('H') )
    h_md    = float( parameters.get('H_md') )
    x_imp   = float( parameters.get('x_imp') )
    rho_ref = float( parameters.get('rho_ref') )
    u_ref   = float( parameters.get('u_ref') )
    p_ref   = float( parameters.get('p_ref') )
    
    dyn_p   = 0.5*rho_ref*u_ref*u_ref
    
    S.df_fric = shift_coordinates( S.df_fric, delta, h_ridge, h_md, x_imp)

    xx     = np.array( S.df_fric['xs'] )
    zz     = np.array( S.df_fric['zs'] )
    fric   = np.array( S.df_fric['fric'] )
    p      = np.array( S.df_pres['p'] )
    p_fluc = np.array( S.df_pres['p`'] )

    npx = len( np.unique(xx) )
    npz = len( np.unique(zz) )

    xx     = xx.reshape( npz, npx )
    zz     = zz.reshape( npz, npx )
    fric   = fric.reshape( npz, npx )
    p      = p.reshape( npz, npx )
    p_fluc = p_fluc.reshape( npz, npx )

    save_isolines( xx, zz, fric,
                   1.0,"projected_separation_line.pkl")
    
    plot_wall_projection( xx, zz, fric/dyn_p*1000.0, filename='fric' )
    
    plot_wall_projection( xx, zz, p/p_ref, filename='pressure' )
    
    plot_wall_projection( xx, zz, p_fluc/p_ref, separation=False, filename='pressure_fluc' )
    
with timer("save spanwise averaged variable distribution along x"):
    
    fric_mean   = np.mean( fric/dyn_p*1000.0, axis=0 )
    p_mean      = np.mean( p/p_ref, axis=0 )
    p_fluc_mean = np.mean( p_fluc/p_ref, axis=0 )

    df_streamwise = pd.DataFrame(columns=['x','fric_ave','p_ave','p_fluc_ave'])
    df_streamwise['x']  = np.unique( xx )
    df_streamwise['Cf'] = np.array( fric_mean )
    df_streamwise['Cp'] = np.array( p_mean )
    df_streamwise['p_fluc'] = np.array( p_fluc_mean )
    
    with open('streamwise_vars.pkl','wb') as f:
        
        pickle.dump( df_streamwise, f )
        
    
         