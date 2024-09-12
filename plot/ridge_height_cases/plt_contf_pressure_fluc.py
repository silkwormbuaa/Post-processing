#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_contf_pressure_fluc.py
@Time    :   2024/03/29 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plot the contour of pressure fluctuation of different ridge height cases
'''

import os
import sys
import pickle
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker
from   scipy.interpolate import griddata

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.plane_analy import compute_DS

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

# =============================================================================
# setting parameters
# =============================================================================

output_dir = '/media/wencan/Expansion/temp/DataPost/lowRe_ridge_height/xy_plane/'
casename   = ['smooth_adiabatic','240211','220927','240210']

df_dirs  = [f'/media/wencan/Expansion/temp/{case}/postprocess/statistics/xy_planes' for case in casename]

label    = ['smooth',  r'$H/\delta_0=0.05$', r'$H/\delta_0=0.10$', r'$H/\delta_0=0.20$',]

varname     = 'p`'         # 'u', 'mach', 'grad_rho', 'T', 'tke', 'p`', 'grad_p', 'DS'
varnorm     =  45447.289   # 45447.289, 1.0
col_map     = 'coolwarm'     #'RdBu_r','Greys_r'
cbar_label  = r'$\sqrt{\langle p ^\prime p ^\prime \rangle}/p_{\infty}$'      # ,'DS' 
cbar_levels = np.linspace(0,0.5,51)  # np.linspace( 0.0, 0.8,33)  np.linspace(0,0.5,51) 
cbar_ticks  = np.linspace(0,0.5,6)   # np.linspace( 0.0, 0.8,5)   np.linspace(0,0.5,6)  
figname0    = 'contour_pressure_fluctuation_awall_linestyle' # 'contour_pressure_fluctuation_awall','contour_DS_awall'
format      = '.pdf'

# =============================================================================

for i in range(2):

    if   i == 0: 
        figname = figname0 + '_ridge'
        y = np.linspace(0.01, 8, 201)
        
    elif i == 1: 
        figname = figname0 + '_valley'
        if casename == 'smooth': y = np.linspace(0.01, 8, 201)
        else:                    y = np.linspace(-0.1, 8, 405)
    
    # create figure
    
    fig, axs = plt.subplots(2,2,figsize=(24,15))
    
    for j, ax in enumerate(axs.flat):
    
        # read in data and interpolate
        os.chdir(df_dirs[j])
        df_slice_file  = f'df_slice_{i:02d}.pkl'
        soniclinefile  = f'soniclines_{i:02d}.pkl'
        seplinefile    = f'seplines_{i:02d}.pkl'
        shockshapefile = f'shockshape_{i:02d}.pkl'
        shockDSfile    = f'shockDS_{i:02d}.pkl'
        
        if not os.path.exists(df_slice_file):
            print(f'\033[93m{casename[j]} {df_slice_file} does not exist!\033[0m')
            continue
        
        df_slice = pd.read_pickle(df_slice_file)
        
        print(f"{casename[j]} is doing interpolation...")
        
        x_slice = np.array( df_slice['xs'] )
        y_slice = np.array( df_slice['y_scale'] )
        
        x = np.linspace(-20,10,301)
        xx,yy = np.meshgrid(x, y)
        
        if varname == 'DS': 
            var_slice = compute_DS(np.array( df_slice['grad_rho'] ), min=0.0, max=2.0 )
            
        else: var_slice  = np.array( df_slice[varname] )
        
        var = griddata( (x_slice,y_slice), 
                        var_slice/varnorm,
                        (xx,yy), method='linear') 
        
        cs = ax.contourf( xx,yy,var,
                          cmap=col_map,
                          levels=cbar_levels,
                          extend='both')
        
        with open(soniclinefile,'rb') as f: lines = pickle.load(f)
        for line in lines:
            ax.plot(line[:,0],line[:,1],'black',linewidth=2.0,linestyle='--')
#            ax.plot(line[:,0],line[:,1],'magenta',linewidth=2.0)
        
        with open(seplinefile,'rb') as f: lines = pickle.load(f)
        for line in lines:
            ax.plot(line[:,0],line[:,1],'black',linewidth=2.5,linestyle=':')
#           ax.plot(line[:,0],line[:,1],'red',linewidth=2.5)
               
        with open(shockshapefile,'rb') as f: lines = pickle.load(f)
        for line in lines:
            ax.plot(line[:,0],line[:,1],'black',linewidth=2.0)
            
        ax.set_xlim([-15,10])
        ax.set_ylim([0,8])
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(2)
        ax.tick_params(axis='both', length=10, width=2, pad=10)
        ax.set_aspect('equal',adjustable='box')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        
        ax.text(1.5,4,label[j])
        
        if j//2 == 1:
            ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
        else: 
            ax.set_xticklabels([])
            
        if j % 2 == 0:
            ax.set_ylabel(r'$y/\delta_0$')
        else:
            ax.set_yticklabels([])
    
    cbar_ax = fig.add_axes([0.4, 0.16, 0.3, 0.03])
    cbar = plt.colorbar(cs,
                        cax=cbar_ax, 
                        orientation='horizontal',
                        ticks=cbar_ticks)
    
    cbar.ax.tick_params( axis='x',
                         direction='in',
                         bottom = True, top=False,
                         length=15.0,
                         width=1.5)

    cbar.ax.set_ylabel(cbar_label, rotation='horizontal', labelpad=0.2 )
    cbar.ax.yaxis.set_label_coords(-0.3,-0.15)
    
    plt.subplots_adjust(wspace=0.10, hspace=0.03, left=0.1, bottom=0.3, right=0.95)

    plt.savefig(output_dir + figname + format)