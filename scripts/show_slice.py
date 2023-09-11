#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_slice.py
@Time    :   2023/09/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   get a slice from statistics.bin
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

from   vista.plane_analy import save_sonic_line

from   vista.plane_analy import save_separation_line

datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'

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
        
        S.compute_vars( block_list, ['mach','RS'] )
        
        S.compute_gradients(block_list, ['schlieren','shadowgraph'],G)
        
with timer("Get slice dataframe "):
    
    df_slice = S.get_slice( G, block_list, indx_slic, 'Z' )

with timer("Interpolate and plot "):
    
    x_slice = np.array( df_slice['x'] )
    y_slice = np.array( df_slice['y'] )
    u_slice = np.array( df_slice['u'] )
    mach_slice = np.array( df_slice['mach'] )
    gradrho_slice = np.array( df_slice['grad_rho'])
    
    x = np.linspace(-50.0,100,301)
    y = np.linspace(0.0, 50.0,101)
    
    xx,yy = np.meshgrid(x,y)
    
#    v = griddata( (x_slice,y_slice), u_slice,
#                  (xx,yy), 
#                  method='linear')
    
    mach = griddata( (x_slice,y_slice), mach_slice,
                     (xx,yy), method='linear')
    
    u    = griddata( (x_slice,y_slice), u_slice,
                     (xx,yy), method='linear')
    
    gradrho = griddata( (x_slice,y_slice), gradrho_slice,
                        (xx,yy), method='linear')
    
    
    save_sonic_line( xx,yy, mach )
    
    save_separation_line( xx, yy, u )
    
    fig,ax = plt.subplots()
    
    ax.contourf( xx,yy, gradrho, levels=51, cmap='bwr')
    
    plt.show()
    
    


