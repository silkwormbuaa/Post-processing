#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_slice.py
@Time    :   2023/09/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   get a slice from statistics.bin (only applicable to 3D)
'''

import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

import numpy             as     np

from   scipy.interpolate import griddata

from   vista.statistic   import StatisticData

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.plane_analy import save_sonic_line

from   vista.plane_analy import save_separation_line

from   vista.plane_analy import shift_coordinates

from   vista.plot_style  import plot_slicez_stat

from   vista.tools       import read_case_parameter

datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'
parametersfile = datapath.split('/results')[0] + '/case_parameters'

# - read in grid info

G = GridData( gridfile )

G.read_grid()

block_list, indx_slic = G.select_sliced_blockgrids( 'Z', 0.0 )

print(f"Selected {len(block_list)} blocks.\n")

# - read statistics data file

S = StatisticData( datafile )

with timer("read selected blocks "):
    
    with open(datafile,'br') as f:
        
        S.read_stat_header( f )
        
        vars = ['u','v','w','T','rho']
        
        S.read_stat_body( f, block_list, vars )
        
        S.compute_vars( block_list, ['mach'] )
        
        S.compute_gradients( block_list, ['schlieren','shadowgraph'],G)
        
with timer("Get slice dataframe "):
    
    df_slice = S.get_slice_df( block_list, G, indx_slic, 'Z' )

with timer("Interpolate and plot "):
    
    parameters = read_case_parameter( parametersfile )
    delta   = float( parameters.get('delta_0') )
    h_ridge = float( parameters.get('H') )
    h_md    = float( parameters.get('H_md') )
    x_imp   = float( parameters.get('x_imp') )
    
    df_slice = shift_coordinates( df_slice, delta, h_ridge, h_md, x_imp )  
    
    x_slice = np.array( df_slice['xs'] )
    y_slice = np.array( df_slice['y_scale'] )
    
    u_slice = np.array( df_slice['u'] )
    mach_slice = np.array( df_slice['mach'] )
    gradrho_slice = np.array( df_slice['grad_rho'] )
    T_slice = np.array( df_slice['T'] )
    
    x = np.linspace(-18,12,301)
    y = np.linspace(0.0,10,101)
    
    xx,yy = np.meshgrid(x,y)
    
    mach = griddata( (x_slice,y_slice), mach_slice,
                     (xx,yy), method='linear')
    
    u    = griddata( (x_slice,y_slice), u_slice,
                     (xx,yy), method='linear')
    
    gradrho = griddata( (x_slice,y_slice), gradrho_slice,
                        (xx,yy), method='linear')

    T    = griddata( (x_slice,y_slice), T_slice,
                     (xx,yy), method='linear')

    
    save_sonic_line( xx,yy, mach )
    
    save_separation_line( xx, yy, u )
    
    cbar = r'$Mach$'
    
    plot_slicez_stat( xx,yy,mach, 
                      filename='MachZ',
                      cbar_label=cbar)

    cbar = r'$<T>/T_{\infty}$'
    
    plot_slicez_stat( xx,yy,T/160.15, 
                      filename='TemperatureZ',
                      cbar_label=cbar)

    cbar = r'$\nabla{\rho}$'
    
    plot_slicez_stat( xx,yy,gradrho, 
                      filename='grad_rho',
                      col_map='Greys',
                      cbar_label=cbar)
    


