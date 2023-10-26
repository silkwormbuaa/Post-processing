# -*- coding: utf-8 -*-
'''
@File    :   pli_p_prime.py
@Time    :   2022/10/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plot spanwise averaged variable distribution along streamwise
'''

import os
import sys
import pickle
from   plt_tools         import PlotDataframe
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker

source_dir = os.path.realpath(__file__).split('plotlines')[0]
sys.path.append( source_dir )

from   vista.line        import LineData

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

OutPath  = "/home/wencanwu/my_simulation/temp/DataPost/streamwise_lines"

sw_pfluc_file = '/home/wencanwu/my_simulation/temp/smooth_wall/line_p_prime.dat'
sw_Cf_file = '/home/wencanwu/my_simulation/temp/smooth_wall/x_cf_STBLI_Wencan.dat'
sw_Cp_file = '/home/wencanwu/my_simulation/temp/smooth_wall/Cf_flat_new.dat'

data1 = '/media/wencanwu/Seagate Expansion Drive/temp/221014/results/wall_projection/streamwise_vars.pkl'
data2 = '/media/wencanwu/Seagate Expansion Drive/temp/220926/results/wall_projection/streamwise_vars.pkl'
data3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/results/wall_projection/streamwise_vars.pkl'
data4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/wall_projection/streamwise_vars.pkl'
data5 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results/wall_projection/streamwise_vars.pkl'


datalist = [data1,   data2,   data3,                   data4,        data5]
dy       = [0.494,   0.468,   0.416,                   0.312,        0.26]
color    = ['green', 'blue', 'black',                  'red',        'purple']
label    = ['2.0',   '1.0',  '0.5',                    '0.25',       '0.125']
lstyle   = [':',     '-.',    (0, (3, 1, 1, 1, 1, 1)), (0, (10, 3)), '-']
width    = [4.0,      4.0,    4.0,                     4.0,          4.0]
lines = []

plt_pwfluc = True
plt_pw     = True
plt_Cf     = True

pure = False

fmt =  '.pdf' # or '.png'

# - read in data files

for i, datafile in enumerate( datalist ):
    
    line = LineData()
    
    with open(datafile,'rb') as f:  line.df = pickle.load( f )
    line.color  = color[i]
    line.label  = r'$\mathrm{D/\delta_0=}$' + label[i]
    line.lstyle = lstyle[i]
    line.width  = width[i]
    
    lines.append( line )

d0_pwfluc = PlotDataframe(sw_pfluc_file)
d0_pw     = PlotDataframe(sw_Cp_file)
d0_f      = PlotDataframe(sw_Cf_file)
d0_f.shift_x( 50.4, 5.2 )
d0_pw.shift_x( 50.4, 5.2 )

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

    ax.plot( d0_pwfluc.df['(x-x_imp)/δ'], 
            d0_pwfluc.df['<p`>_'],
            'gray', 
            label=r'$smooth$', 
            ls   ='--',
            linewidth=4)

    # separation locations

    #ax.plot( [-10.6786859,-10.6786859],
    #         [0.0,0.1],
    #         'purple',
    #         linewidth=1)
    #
    #ax.plot( [-10.6574429,-10.6574429],
    #         [0.0,0.1],
    #            'red',
    #            linewidth=1)
    #
    #ax.plot( [-9.179693795,-9.179693795],
    #         [0.0,0.1],
    #            'black',
    #            linewidth=1)
    #
    #ax.plot( [-8.316817364,-8.316817364],
    #         [0.0,0.1],
    #            'blue',
    #            linewidth=1)
    #
    #ax.plot( [-8.405281852,-8.405281852],
    #         [0.0,0.1],
    #            'green',
    #            linewidth=1)
    #
    #ax.plot( [-8.56077913,-8.56077913],
    #         [0.0,0.1],
    #            'gray',
    #            linewidth=1)


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
    ax.set_ylim([0.02,0.1])
    ax.grid(visible=True, which='both',axis='both',color='gray',
                linestyle='--',linewidth=0.2)

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

        ax.set_ylabel(r"$\sqrt{<p'p'>}/p_{\infty}$" )
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

    ax.plot( d0_pw.df['x_s'], 
            d0_pw.df['Cp'],
            'gray', 
            label=r'$smooth$', 
            ls   ='--',
            linewidth=4)

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
    ax.grid(visible=True, which='both',axis='both',color='gray',
                linestyle='--',linewidth=0.2)

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

    ax.plot( d0_f.df['x_s'], 
            d0_f.df['Cf']*1000,
            'gray', 
            label=r'$smooth$', 
            ls   ='--',
            linewidth=4)
    
    ax.plot( [-20,12],
             [0,0],
             'black',
             ls = '--',
             linewidth=1 )

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
    ax.set_ylim([-2.5,4.5])
    
    ax.grid(visible=True, which='both',axis='both',color='gray',
                linestyle='--',linewidth=0.2)

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
        
        ax.set_ylabel("$C_fx10^3$")
        ax.tick_params(axis='y', pad=10)

#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
        
    plt.savefig( figname + fmt )
    plt.show()