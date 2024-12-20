#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_streamwise_vars.py
@Time    :   2024/08/16 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   plotting streamwise variable comparison for different Re cases
'''

import os
import sys
import pickle
import pandas            as     pd
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import LineData

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

OutPath  = "/media/wencan/Expansion/temp/DataPost/midRe/averaged_streamwise_vars_new"

data0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data1 = '/media/wencan/Expansion/temp/smooth_mid/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data2 = '/media/wencan/Expansion/temp/220927/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data3 = '/media/wencan/Expansion/temp/231124/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data4 = '/media/wencan/Expansion/temp/241018/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data5 = '/media/wencan/Expansion/temp/241030/postprocess/statistics/wall_projection/streamwise_vars.pkl'

sw_iw_pw_file = '/media/wencan/Expansion/temp/data_luis/midRe_streamwise/pw_SBLI_B2.dat'
sw_iw_pwrms_file = '/media/wencan/Expansion/temp/data_luis/midRe_streamwise/pw_rms_SBLI_B2.dat'
sw_iw_Cf_file = '/media/wencan/Expansion/temp/data_luis/midRe_streamwise/Cf_SBLI_B2.dat'

add_sw_luis = False

datalist = [data0,   data1,   data2,   data3, data4, data5]
color    = ['gray', 'blue',   'gray', 'pink', 'red', 'green'] 
label    = ['lowRe_smooth', 'midRe_smooth', 'lowRe_rough', 'midRe_0.1', 'midRe_0.2','midRe_0.026']
lstyle   = ['--',   '--',    '-.',  '-.', '-.','-.']
width    = [4.0,     4.0,    4.0,    4.0, 4.0, 4.0]
lines = []

plt_pwfluc = True
plt_pw     = True
plt_Cf     = True

figsizes = [15,8]

pure = False
show_label = True

fmt =  'sw.png' # or '.png'

# - read in data files

for i, datafile in enumerate( datalist ):
    
    line = LineData()
    
    with open(datafile,'rb') as f:  line.df = pickle.load( f )
    line.color  = color[i]
    line.label  = label[i]
    line.lstyle = lstyle[i]
    line.width  = width[i]
    
    lines.append( line )

os.chdir(OutPath)

# ----------------------------------------------------------------------
# >>> plot wall pressure fluctuation                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_pwfluc:

    fig, ax = plt.subplots( figsize=figsizes, constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                 line.df['p_fluc'],
                 line.color,
                 ls = line.lstyle,
                 label = line.label,
                 linewidth = line.width)
    
    if add_sw_luis:
        swdf = pd.read_csv(sw_iw_pwrms_file, delimiter=r'\s+')
        ax.plot( swdf['x']*7.15/5.2, 
                 swdf['y'],
                 'black',
                 ls = '-',
                 label = 'midRe_smooth',
                 linewidth = 4.0)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=15,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1)

    ax.set_xlim([-20,10])
    ax.set_ylim([0.01,0.12])
#    ax.grid(visible=True, which='both',axis='both',color='gray',
#                linestyle='--',linewidth=0.2)

    figname = "pressure_fluctuation"

    # Adjust the spacing around the plot to remove the white margin
    if pure:
        figname += '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1) 
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        
        ax.set_xlabel(r"$(x-x_{imp})/\delta_0$", labelpad=-5 )  
        ax.tick_params(axis='x', pad=15)

        ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_{\infty}$" )
        ax.tick_params(axis='y', pad=10)
        
        if show_label:
            ax.legend( fontsize=20 ) 
        
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()
    
# ----------------------------------------------------------------------
# >>> plot wall pressure                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_pw:

    fig, ax = plt.subplots( figsize=figsizes, constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                line.df['Cp'],
                line.color,
                ls = line.lstyle,
                label = line.label,
                linewidth = line.width)

    if add_sw_luis:
        swdf = pd.read_csv(sw_iw_pw_file, delimiter=r'\s+')
        ax.plot( swdf['x']*7.15/5.2, 
                 swdf['y'],
                 'black',
                 ls = '-',
                 label = 'midRe_smooth',
                 linewidth = 4.0)
        
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=15,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1)

    ax.set_xlim([-20.0,10.0])
    ax.set_ylim([0.8,2.5])
#    ax.grid(visible=True, which='both',axis='both',color='gray',
#                linestyle='--',linewidth=0.2)

    figname = "wall_pressure"

    # Adjust the spacing around the plot to remove the white margin
    if pure:
        figname += '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1) 
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        
        ax.set_xlabel("$(x-x_{imp})/\delta_0$", labelpad=-5 )
        ax.tick_params(axis='x', pad=15)
        
        ax.set_ylabel("$<p_w>/p_{\infty}$")
        ax.tick_params(axis='y', pad=10)

        if show_label:
            ax.legend( fontsize=30 ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)        
        
    plt.savefig( figname + fmt )
    plt.show()
    
# ----------------------------------------------------------------------
# >>> plot friction coefficient                                   (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Cf:
    
    fig, ax = plt.subplots( figsize=figsizes, constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                line.df['Cf'],
                line.color,
                ls = line.lstyle,
                label = line.label,
                linewidth = line.width)

    if add_sw_luis:
        swdf = pd.read_csv(sw_iw_Cf_file, delimiter=r'\s+')
        ax.plot( swdf['x']*7.15/5.2, 
                 swdf['y']*1000.0,
                 'black',
                 ls = '-',
                 label = 'midRe_smooth',
                 linewidth = 4.0)

    ax.plot( [-20,12],
             [0,0],
             'black',
             ls = '--',
             linewidth=2 )
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=15,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1)

    ax.set_xlim([-20.0,10.0])
    ax.set_ylim([-2.5,5.0])
    
#    ax.grid(visible=True, which='both',axis='both',color='gray',
#                linestyle='--',linewidth=0.2)

    figname = "friction_coefficient"

    # Adjust the spacing around the plot to remove the white margin
    if pure:
        figname += '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1) 
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        
        ax.set_xlabel("$(x-x_{imp})/\delta_0$", labelpad=-5) 
        ax.tick_params(axis='x', pad=15)
        
        ax.set_ylabel(r"$C_f \times 10^3$")
        ax.tick_params(axis='y', pad=10)

        if show_label:
            ax.legend( fontsize=30 ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
        
    plt.savefig( figname + fmt )
    plt.show()