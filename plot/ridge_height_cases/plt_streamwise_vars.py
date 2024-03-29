#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_streamwise_vars_H.py
@Time    :   2024/03/05 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   plotting streamwise variables for different ridge heights cases
'''

import os
import sys
import pickle
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import LineData
from   vista.plot_tools  import PlotDataframe

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

OutPath  = "/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe_ridge_height/averaged_streamwise_vars"

# sw_pfluc_file = '/home/wencanwu/my_simulation/temp/smooth_wall/line_p_prime.dat'
# sw_Cf_file = '/home/wencanwu/my_simulation/temp/smooth_wall/x_cf_STBLI_Wencan.dat'
# sw_Cp_file = '/home/wencanwu/my_simulation/temp/smooth_wall/Cf_flat_new.dat'

data1 = '/media/wencanwu/Seagate Expansion Drive1/temp/240211/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data2 = '/media/wencanwu/Seagate Expansion Drive1/temp/220927/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data3 = '/media/wencanwu/Seagate Expansion Drive1/temp/240210/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data4 = '/media/wencanwu/Seagate Expansion Drive1/temp/smooth_adiabatic/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data5 = '/media/wencanwu/Seagate Expansion Drive1/temp/smooth_isothermal/postprocess/statistics/wall_projection/streamwise_vars.pkl'

datalist = [data1,   data2,   data3, data4]   #, data5 ]
color    = ['black',  'red', 'blue','gray'] #, 'green' ]
label    = ['0.05',   '0.1',  '0.2', 'smooth_awall','smooth_iwall']
lstyle   = [':',     '-.',    (0, (3, 1, 1, 1, 1, 1)), 'dashed', 'dotted']
width    = [4.0,      4.0,    4.0 , 4.0, 4.0 ]
lines = []

plt_pwfluc = True
plt_pw     = True
plt_Cf     = True

pure = False

fmt =  '.png' # or '.png'

# - read in data files

for i, datafile in enumerate( datalist ):
    
    line = LineData()
    
    with open(datafile,'rb') as f:  line.df = pickle.load( f )
    line.color  = color[i]
    line.label  = r'$\mathrm{H/\delta_0=}$' + label[i]
    line.lstyle = lstyle[i]
    line.width  = width[i]
    
    lines.append( line )

# d0_pwfluc = PlotDataframe(sw_pfluc_file)
# d0_pw     = PlotDataframe(sw_Cp_file)
# d0_f      = PlotDataframe(sw_Cf_file)
# d0_f.shift_x( 50.4, 5.2 )
# d0_pw.shift_x( 50.4, 5.2 )

os.chdir(OutPath)

# ----------------------------------------------------------------------
# >>> plot wall pressure fluctuation                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/05  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_pwfluc:

    fig, ax = plt.subplots( figsize=[15,5], constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                 line.df['p_fluc'],
                 line.color,
                 ls = line.lstyle,
                 label = line.label,
                 linewidth = line.width)

    # ax.plot( d0_pwfluc.df['(x-x_imp)/δ'], 
    #         d0_pwfluc.df['<p`>_'],
    #         'gray', 
    #         label=r'$smooth$', 
    #         ls   ='--',
    #         linewidth=4)

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
    ax.set_ylim([0.01,0.09])
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

#        ax.legend( ) 
        
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
# 2023/10/05  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_pw:

    fig, ax = plt.subplots( figsize=[15,5], constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                line.df['Cp'],
                line.color,
                ls = line.lstyle,
                label = line.label,
                linewidth = line.width)

    # ax.plot( d0_pw.df['x_s'], 
    #         d0_pw.df['Cp'],
    #         'gray', 
    #         label=r'$smooth$', 
    #         ls   ='--',
    #         linewidth=4)

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

#        ax.legend( ) 

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
# 2023/10/05  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Cf:
    
    fig, ax = plt.subplots( figsize=[15,5], constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                line.df['Cf'],
                line.color,
                ls = line.lstyle,
                label = line.label,
                linewidth = line.width)

    # ax.plot( d0_f.df['x_s'], 
    #         d0_f.df['Cf']*1000,
    #         'gray', 
    #         label=r'$smooth$', 
    #         ls   ='--',
    #         linewidth=4)
    
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

#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
        
    plt.savefig( figname + fmt )
    plt.show()