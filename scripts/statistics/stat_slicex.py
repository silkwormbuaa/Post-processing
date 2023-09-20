#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_slicex.py
@Time    :   2023/09/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

import numpy             as     np

import matplotlib.pyplot as     plt

from   scipy.interpolate import griddata

from   vista.statistic   import StatisticData

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.tools       import define_wall_shape

from   vista.tools       import read_case_parameter

from   vista.plane_analy import save_sonic_line

from   vista.plane_analy import shift_coordinates

from   vista.plot_style  import plot_slicex_stat

loc = -57.0


datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'
parametersfile = datapath.split('/results')[0] + '/case_parameters'

# - read in grid info

G = GridData( gridfile )

G.read_grid()

block_list, indx_slic = G.select_sliced_blockgrids( 'X', loc )

print(f"Selected {len(block_list)} blocks.\n")

# - read statistics data file

S = StatisticData( datafile )

with timer("read selected blocks "):
    
    with open(datafile,'br') as f:
        
        S.read_stat_header( f )
        
        vars = ['u','v','w','uu','vv','ww','uv','vw','T']
        
        S.read_stat_body( f, block_list, vars )
        
        S.compute_vars( block_list, ['mach','RS'] )
        
        S.compute_gradients( block_list, ['vorticity'],G,)
        
        S.compute_source_terms( block_list, G )
        
with timer("Get slice dataframe and match grids"):
    
    df_slice = S.get_slice_df( block_list, G, indx_slic, 'X' )
    
with timer("Interpolate and plot "):
    
    parameters = read_case_parameter( parametersfile )
    delta   = float( parameters.get('delta_0') )
    h_ridge = float( parameters.get('H') )
    h_md    = float( parameters.get('H_md') )
    x_imp   = float( parameters.get('x_imp') )
    casecode =  str( parameters.get('casecode') )
    
    df_slice = shift_coordinates( df_slice, delta, h_ridge, h_md, x_imp )  
    
    y_slice = np.array( df_slice['ys'] )
    z_slice = np.array( df_slice['zs'] )
    
    mach_slice = np.array( df_slice['mach'] )
    tke_slice = np.array( df_slice['tke'] )
    S_slice = np.array( df_slice['S'] )
    w1_slice = np.array( df_slice['w1'] )
    u_slice = np.array( df_slice['u'] )
    v_slice = np.array( df_slice['v'] )
    w_slice = np.array( df_slice['w'] )
    
    z = np.linspace(-1.0,1.0, 201)
    y = np.linspace(-0.1, 1.1, 121)
    
    zz,yy = np.meshgrid(z,y)
    
    mach = griddata( (z_slice,y_slice), mach_slice,
                     (zz,yy), method='linear')
    
    tke = griddata( (z_slice,y_slice), tke_slice,
                    (zz,yy), method='linear')
    
    S = griddata( (z_slice,y_slice), S_slice,
                   (zz,yy), method='linear')
    
    w1 = griddata( (z_slice,y_slice), w1_slice,
                   (zz,yy), method='linear')    

    u = griddata( (z_slice,y_slice), u_slice,
                   (zz,yy), method='linear') 

    v = griddata( (z_slice,y_slice), v_slice,
                   (zz,yy), method='linear') 

    w = griddata( (z_slice,y_slice), w_slice,
                   (zz,yy), method='linear')   
    
    
    fig, ax = plt.subplots()
    
#    ax.contourf( zz, yy, mach, levels=51, cmap='bwr' )

    save_sonic_line( zz, yy, mach )
    
    define_wall_shape( z*5.2, casecode=casecode, yshift=(h_ridge-h_md) )
    
    cbar = r'$Mach$'
    plot_slicex_stat( zz, yy, mach,
                      filename='MachX',
                      cbar_label=cbar,
                      col_map='coolwarm')

    cbar = r'$S$'
    cbar_levels=np.linspace(-20000.0,20000.,51)
    plot_slicex_stat( zz, yy, S,
                      filename='S_X',
                      cbar_label=cbar,
                      cbar_levels=cbar_levels)
    
    cbar = r'$\omega_x$'
    plot_slicex_stat( zz, yy, w1,
                      filename='vorticity_X',
                      cbar_label=cbar,
                      col_map='coolwarm')

    cbar = r'$u$'
    cbar_levels=np.linspace(0.0,507.0,51)
    plot_slicex_stat( zz, yy, u,
                      filename='u',
                      cbar_label=cbar,
                      col_map='coolwarm',
                      cbar_levels=cbar_levels)

    cbar = r'$v$'
    plot_slicex_stat( zz, yy, v,
                      vectors=[w,v],
                      filename='vertical_velocity_X',
                      cbar_label=cbar,
                      col_map='coolwarm')
    
    cbar = r'$tke$'
    plot_slicex_stat( zz, yy, tke,
                      vectors=[w,v],
                      filename='tke_X',
                      cbar_label=cbar,
                      col_map='coolwarm')