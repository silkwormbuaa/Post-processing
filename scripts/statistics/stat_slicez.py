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

import numpy             as     np

from   scipy.interpolate import griddata

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.statistic   import StatisticData

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.plane_analy import save_sonic_line
from   vista.plane_analy import save_separation_line
from   vista.plane_analy import shift_coordinates

from   vista.plot_style  import plot_slicez_stat

from   vista.tools       import read_case_parameter



# =============================================================================

# locs = [ridge location, valley location]

outfoler = '/xy_planes'

# =============================================================================


datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'
outpath  = datapath + outfoler
parametersfile = datapath.split('/results')[0] + '/case_parameters'

# - read in case parameters

parameters = read_case_parameter( parametersfile )
delta   = float( parameters.get('delta_0') )
h_ridge = float( parameters.get('H') )
h_md    = float( parameters.get('H_md') )
x_imp   = float( parameters.get('x_imp') )
p_ref   = float( parameters.get('p_ref') )
u_ref   = float( parameters.get('u_ref') )
casecode =  str( parameters.get('casecode') )


locs = [ 0.0, 0.5*float( parameters.get('D')) ]

# - read in grid info

G = GridData( gridfile )

G.read_grid()

# - enter outpath

if not os.path.exists(outpath): 
    os.mkdir( outpath )
    print(f"Created directory {outpath}.\n")

os.chdir(outpath)

# - do slicing and output slices

for i, loc in enumerate( locs ):

    block_list, indx_slic = G.select_sliced_blockgrids( 'Z', loc )

    print(f"Selected {len(block_list)} blocks.\n")

    # - read statistics data file

    S = StatisticData( datafile )

    with timer("read selected blocks "):
        
        with open(datafile,'br') as f:
            
            S.read_stat_header( f )
            
            vars = ['u','v','w','T','rho','uu','vv','ww','uv','pp','p']
            
            S.read_stat_body( f, block_list, vars )
            
            S.compute_vars( block_list, ['mach','RS','p`'] )
            
            S.compute_gradients( block_list, 
                                ['schlieren','shadowgraph','vorticity','grad_p'],
                                G)
            
    with timer("Get slice dataframe "):
        
        df_slice = S.get_slice_df( block_list, G, indx_slic, 'Z' )

    with timer("Interpolate and plot "):
        
        df_slice = shift_coordinates( df_slice, delta, h_ridge, h_md, x_imp )  
        
        x_slice = np.array( df_slice['xs'] )
        y_slice = np.array( df_slice['y_scale'] )
        
        u_slice = np.array( df_slice['u'] )
        mach_slice = np.array( df_slice['mach'] )
        gradrho_slice = np.array( df_slice['grad_rho'] )
        T_slice = np.array( df_slice['T'] )
        tke_slice = np.array( df_slice['tke'] )
        p_fluc_slice = np.array( df_slice['p`'] )
        grad_p_slice = np.array( df_slice['grad_p'] )
        
        x = np.linspace(-20,10,301)
        
        if i == 0:
            y = np.linspace(0.0, 8, 201)
            wall = False
        
        if i == 1:
            if casecode == 'smooth_wall':  y = np.linspace(0.0, 8, 201)
            else:                          y = np.linspace(-0.1,8, 405)
            wall=False
        
        xx,yy = np.meshgrid(x,y)
        
        mach = griddata( (x_slice,y_slice), mach_slice,
                        (xx,yy), method='linear')
        
        u    = griddata( (x_slice,y_slice), u_slice,
                        (xx,yy), method='linear')
        
        gradrho = griddata( (x_slice,y_slice), gradrho_slice,
                            (xx,yy), method='linear')

        T    = griddata( (x_slice,y_slice), T_slice,
                        (xx,yy), method='linear')

        tke   = griddata( (x_slice,y_slice), tke_slice,
                        (xx,yy), method='linear')

        p_fluc = griddata( (x_slice,y_slice), p_fluc_slice,
                        (xx,yy), method='linear') 

        grad_p = griddata( (x_slice,y_slice), grad_p_slice,
                        (xx,yy), method='linear') 
        
        save_sonic_line( xx,yy, mach )
        
        save_separation_line( xx, yy, u )
        
        cbar = r'$Mach$'
        
        plot_slicez_stat( xx,yy,mach, 
                          filename='MachZ_'+str(i+1),
                          cbar_label=cbar,
                          col_map='coolwarm',
                          wall=wall)

        cbar = r'$<T>/T_{\infty}$'
        
        plot_slicez_stat( xx,yy,T/160.15, 
                          filename='TemperatureZ_'+str(i+1),
                          col_map='plasma',
                          cbar_label=cbar,
                          wall=wall)

        cbar = r'$\nabla{\rho}$'
        
        plot_slicez_stat( xx,yy,gradrho, 
                        filename='grad_rho_'+str(i+1),
                        col_map='Greys',
                        cbar_label=cbar,
                        wall=wall)

        cbar = r'$tke$'
        
        plot_slicez_stat( xx,yy,tke, 
                        filename='tke_'+str(i+1),
                        col_map='coolwarm',
                        cbar_label=cbar,
                        wall=wall)
        
        cbar = 'pressure fluctuation'
        
        plot_slicez_stat( xx,yy,p_fluc, 
                        filename='pressure_fluc_'+str(i+1),
                        col_map='coolwarm',
                        cbar_label=cbar,
                        wall=wall)

        cbar = 'pressure gradient'
        
        plot_slicez_stat( xx,yy,grad_p*delta/p_ref, 
                        filename='p_grad_'+str(i+1),
                        col_map='coolwarm',
                        cbar_label=cbar,
                        wall=wall)
    


