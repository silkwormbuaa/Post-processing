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
from   vista.plane_analy import periodic_average

from   vista.plot_style  import plot_slicex_stat

# =============================================================================

locs_delta = np.linspace(-20,-20,1)
outfolder  = '/yz_planes'

# =============================================================================

datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'
outpath  = datapath + outfolder
parametersfile = datapath.split('/results')[0] + '/case_parameters'

# - read in case paramters

parameters = read_case_parameter( parametersfile )
delta   = float( parameters.get('delta_0') )
h_ridge = float( parameters.get('H') )
h_md    = float( parameters.get('H_md') )
x_imp   = float( parameters.get('x_imp') )
p_ref   = float( parameters.get('p_ref') )
u_ref   = float( parameters.get('u_ref') )
casecode =  str( parameters.get('casecode') )

locs = locs_delta*delta + 50.4

# - read in grid info

G = GridData( gridfile )

G.read_grid()

# - enter outpath

if not os.path.exists(outpath): 
    os.mkdir( outpath )
    print(f"Created directory {outpath}.\n")

os.chdir(outpath)

# do slicing and output slices

for i, loc in enumerate(locs):

    loc_delta = locs_delta[i]
    title = 'x= '+str(loc_delta)
    
    print(f"Start doing slicing at x = {loc_delta:10.2f}.\n")

    block_list, indx_slic = G.select_sliced_blockgrids( 'X', loc )

    print(f"Selected {len(block_list)} blocks.\n")

    # - read statistics data file

    S = StatisticData( datafile )

    with timer("read selected blocks "):
        
        with open(datafile,'br') as f:
            
            S.read_stat_header( f )
            
            vars = ['u','v','w','uu','vv','ww','uv','vw','T','p','pp']
            
            S.read_stat_body( f, block_list, vars )
            
            S.compute_vars( block_list, ['mach','RS','p`'] )
            
            S.compute_gradients( block_list, ['vorticity'], G)
            
            S.compute_source_terms( block_list, G )
            
    with timer("Get slice dataframe and match grids"):
        
        df_slice = S.get_slice_df( block_list, G, indx_slic, 'X' )
        
    with timer("Interpolate and plot "):
               
        df_slice = shift_coordinates( df_slice, delta, h_ridge, h_md, x_imp )  
        
        y_slice = np.array( df_slice['ys'] )
        z_slice = np.array( df_slice['zs'] )
        
        mach_slice = np.array( df_slice['mach'] )
        tke_slice = np.array( df_slice['tke'] )
        RS_slice = np.array( df_slice['u`v`'] )
        S_slice = np.array( df_slice['S'] )
        w1_slice = np.array( df_slice['w1'] )
        u_slice = np.array( df_slice['u'] )
        v_slice = np.array( df_slice['v'] )
        w_slice = np.array( df_slice['w'] )
        p_fluc_slice = np.array( df_slice['p`'] )
        
        # generate interpolation grid
        
        z = np.linspace(-1.0,1.0, 320)
        if casecode == 'smooth_wall':
            y = np.linspace(0.02, 1.1, 55)
        else:
            y = np.linspace(-0.1, 1.1, 121)
        
        zz,yy = np.meshgrid(z,y)
        
        # mapping variables
        
        mach = griddata( (z_slice,y_slice), mach_slice,
                         (zz,yy), method='linear')
        
        tke = griddata( (z_slice,y_slice), tke_slice,
                        (zz,yy), method='linear')

        RS = griddata( (z_slice,y_slice), RS_slice,
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

        p_fluc = griddata( (z_slice,y_slice), p_fluc_slice,
                           (zz,yy), method='linear')   
        
# ------ extending corner grid for smooth wall

        if casecode == 'smooth_wall':
            
            len_ext = np.shape(zz)[1]
            zz = np.concatenate(([zz[0,:]],zz),axis=0)
            yy = np.concatenate(([np.zeros(len_ext)],yy),axis=0)
            mach = np.concatenate(([np.zeros(len_ext)],mach),axis=0)
            tke = np.concatenate(([np.zeros(len_ext)],tke),axis=0)
            RS = np.concatenate(([np.zeros(len_ext)],RS),axis=0)
            S  = np.concatenate(([np.zeros(len_ext)],S),axis=0)
            w1 = np.concatenate(([w1[0,:]],w1),axis=0)
            u  = np.concatenate(([np.zeros(len_ext)],u),axis=0)
            v  = np.concatenate(([np.zeros(len_ext)],v),axis=0)
            w  = np.concatenate(([np.zeros(len_ext)],w),axis=0)
            p_fluc = np.concatenate(([p_fluc[0,:]],p_fluc),axis=0)
        
        fig, ax = plt.subplots()
        
    #    ax.contourf( zz, yy, mach, levels=51, cmap='bwr' )

        save_sonic_line( zz, yy, mach )
        
        define_wall_shape( z*5.2, casecode=casecode, yshift=(h_ridge-h_md) )
        
        # Do periodic average for smallest ridge spacing case
        
        if casecode == "1221":  # otherwise streamline looks too messy
            w = periodic_average(w,16,axis=1)
            v = periodic_average(v,16,axis=1)
            mach   = periodic_average(mach,16,axis=1)
            tke    = periodic_average(tke,16,axis=1)
            RS     = periodic_average(RS,16,axis=1)
            w1     = periodic_average(w1,16,axis=1)
            p_fluc = periodic_average(p_fluc,16,axis=1)
        
        cbar = r'$\langle Mach \rangle$'
        cbar_levels = np.linspace(0.0,2.0,21)
        cbar_ticks  = np.linspace(0.0,2.0,5)
        if casecode == 'smooth_wall':
            plot_slicex_stat( zz, yy, mach,
                              filename='Mach_'+str(i+1),
                              col_map='coolwarm',
                              cbar_label=cbar,
                              cbar_levels=cbar_levels,
                              cbar_ticks=cbar_ticks,
                              title=title)
        else:
            plot_slicex_stat( zz, yy, mach,
                              vectors=[w,v],
                              filename='Mach_'+str(i+1),
                              col_map='coolwarm',
                              cbar_label=cbar,
                              cbar_levels=cbar_levels,
                              cbar_ticks=cbar_ticks,
                              title=title)

#        cbar = r'$S$'
#        cbar_levels=np.linspace(-20000.0,20000.,51)
#        plot_slicex_stat( zz, yy, S,
#                          filename='S_'+str(i+1),
#                          cbar_label=cbar,
#                          cbar_levels=cbar_levels,
#                          title=title)
        
        cbar = r'$\frac{\langle \omega_x \rangle \delta_0}{u_\infty}$ '
        cbar_levels = np.linspace(-1.0,1.0,21)
        cbar_ticks  = np.linspace(-1.0,1.0,5)
        plot_slicex_stat( zz, yy, w1*delta/u_ref,
                          vectors=[w,v],
                          arrow=True,
                          filename='vorticity_'+str(i+1),
                          col_map='coolwarm',
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          title=title)

        cbar = r'$\langle u \rangle /u_{\infty}$'
        cbar_levels = np.linspace(0.0,1.0,21)
        cbar_ticks  = np.linspace(0.0,1.0,5)
        plot_slicex_stat( zz, yy, u/u_ref,
                          filename='u_'+str(i+1),
                          cbar_label=cbar,
                          col_map='coolwarm',
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          title=title)

        cbar = r'$\langle v \rangle$'
        cbar_levels = np.linspace(-3.0,3.0,31)
        cbar_ticks  = np.linspace(-3.0,3.0,5)
        plot_slicex_stat( zz, yy, v/u_ref*100,
                          vectors=[w,v],
                          arrow=True,
                          filename='v_'+str(i+1),
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          col_map='coolwarm',
                          title=title)
        
        cbar = r'$\langle tke ^+ \rangle \cdot 100$'
        cbar_levels = np.linspace(0.0,2.5,26)
        cbar_ticks  = np.linspace(0.0,2.5,6)
        plot_slicex_stat( zz, yy, tke/(u_ref**2)*100,
                          filename='tke_'+str(i+1),
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          col_map='coolwarm',
                          title=title)

        cbar = r"$-\langle u^{'} v^{'} \rangle$"
        plot_slicex_stat( zz, yy, -RS/(u_ref**2),
                          filename='RS_'+str(i+1),
                          cbar_label=cbar,
                          col_map='coolwarm',
                          title=title)

        cbar = r"$\sqrt{\langle p^{'}p^{'}\rangle}\cdot 100$"
        cbar_levels = np.linspace(0.0, 3.2,33)
        cbar_ticks  = np.linspace(0.0,3.2,5)
        plot_slicex_stat( zz, yy, p_fluc/p_ref*100,
                          filename='p_fluc_'+str(i+1),
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          col_map='coolwarm',
                          cbar_ticks=cbar_ticks,
                          title=title)
        
    print(f"Finished doing slicing at x = {loc_delta:10.2f}.",end='') 
    print(f" Progress {(i+1)/len(locs)*100:5.2f} %.")
    print(f"=====================================\n")