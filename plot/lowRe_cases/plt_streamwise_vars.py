# -*- coding: utf-8 -*-
'''
@File    :   plt_streamwise_vars.py
@Time    :   2022/10/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plot spanwise averaged variable distribution along streamwise direction for low Re cases
'''

import os
import sys
import pickle
import numpy             as     np
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import LineData
from   vista.tools       import create_linear_interpolator


plt.rcParams["text.usetex"]         = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family']         = "Times New Roman"
plt.rcParams['font.size']           = 40

OutPath  = "/media/wencan/Expansion/temp/DataPost/lowRe/averaged_streamwise_vars"

# sw_pfluc_file = '/home/wencanwu/my_simulation/temp/smooth_wall/line_p_prime.dat'
# sw_Cf_file = '/home/wencanwu/my_simulation/temp/smooth_wall/x_cf_STBLI_Wencan.dat'
# sw_Cp_file = '/home/wencanwu/my_simulation/temp/smooth_wall/Cf_flat_new.dat'

data1 = '/media/wencan/Expansion/temp/221014/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data2 = '/media/wencan/Expansion/temp/220926/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data3 = '/media/wencan/Expansion/temp/220825/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data4 = '/media/wencan/Expansion/temp/220927/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data5 = '/media/wencan/Expansion/temp/221221/postprocess/statistics/wall_projection/streamwise_vars.pkl'

data6 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/wall_projection/streamwise_vars.pkl'

datalist = [data1,   data2,   data3,                   data4,        data5,     data6]
dy       = [0.494,   0.468,   0.416,                   0.312,        0.26,      0.0]
color    = ['green', 'blue', 'black',                  'red',        'purple',  'gray']
label    = ['2.0',   '1.0',  '0.5',                    '0.25',       '0.125',   'awall']
lstyle   = [':',     '-.',    (0, (3, 1, 1, 1, 1, 1)), (0, (10, 3)), '-',       '--']
width    = [4.0,      4.0,    4.0,                     4.0,          4.0,       4.0]

x_sep    = [-8.405,  -8.317,  -9.18,                  -10.84,       -10.67,  -8.42 ]
x_att    = [1.400 ,  1.360 ,  2.16 ,                    2.45,       2.59  ,  1.10  ]
x_pfmax  = [-7.712,  -7.49 ,  -8.71,                  -10.06,       -8.93 ,  -7.136]

lines    = []

plt_pwfluc    = False
plt_pwfluc_ln = False
plt_pw        = False
plt_pwg       = True
plt_Cf        = True

pure          = False
show_label    = False

figsize       = [15,8]
fmt           =  '.png' # or '.png'

# - read in data files

for i, datafile in enumerate( datalist ):
    
    line = LineData()
    
    with open(datafile,'rb') as f:  line.df = pickle.load( f )
    line.color  = color[i]
    line.label  = r'$\mathrm{D/\delta_0=}$' + label[i]
    line.lstyle = lstyle[i]
    line.width  = width[i]
    
    lines.append( line )


os.chdir(OutPath)

# adjust plotting style

def adjust_plotting( ax:plt.Axes ):

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
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))

    ax.set_xlabel(r"$(x-x_{imp})/\delta_0$", labelpad=-5 )  
    ax.tick_params(axis='x', pad=15)
    
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    # ax.xaxis.set_ticklabels([])
    # ax.yaxis.set_ticklabels([])

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

    fig, ax = plt.subplots( figsize=figsize, constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                 line.df['p_fluc'],
                 line.color,
                 ls = line.lstyle,
                 label = line.label,
                 linewidth = line.width)

    # to avoid line overlapped on marker
    
    for i,line in enumerate( lines ):
        interpolator = create_linear_interpolator( line.df['x'], line.df['p_fluc'])
        ax.plot( x_sep[i],   interpolator(x_sep[i]), 'p',   color=line.color, ms=15)
        ax.plot( x_att[i],   interpolator(x_att[i]), 'p',   color=line.color, ms=15)
        ax.plot( x_pfmax[i], interpolator(x_pfmax[i]), '*', color=line.color, ms=15)

    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
    ax.set_ylim([0.02,0.1])
    ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_{\infty}$" )
    ax.tick_params(axis='y', pad=10)
    
    adjust_plotting( ax )
    
    figname = "pressure_fluctuation_awall"

    plt.savefig( figname + fmt )
    print( f"{figname.ljust(30)} is output in {OutPath}." )
    plt.show()


# ----------------------------------------------------------------------
# >>> plot locally normalized wall pressure fluctuations         (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/02/05  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_pwfluc_ln:
    
    fig, ax = plt.subplots( figsize=figsize, constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                 line.df['p_fluc']/line.df['Cp'],
                 line.color,
                 ls = line.lstyle,
                 label = line.label,
                 linewidth = line.width)
    
    # to avoid line overlapped on marker
    
    for i,line in enumerate( lines ):
        interpolator = create_linear_interpolator( line.df['x'], line.df['p_fluc']/line.df['Cp'])
        ax.plot( x_sep[i],   interpolator(x_sep[i]),   'p', color=line.color, ms=15)
        ax.plot( x_att[i],   interpolator(x_att[i]),   'p', color=line.color, ms=15)
        ax.plot( x_pfmax[i], interpolator(x_pfmax[i]), '*', color=line.color, ms=15)

    # adjust the y axis

    ax.set_ylim([0.01,0.08])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
    ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_{w}$" )
    ax.tick_params(axis='y', pad=10)
    
    adjust_plotting( ax )
    
    if show_label:
        ax.legend( fontsize=20 ) 
        
    figname = "pressure_fluctuation_ln"
    plt.savefig( figname + fmt )
    print( f"{figname.ljust(30)} is output in {OutPath}." )
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

    fig, ax = plt.subplots( figsize=figsize, constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                line.df['Cp'],
                line.color,
                ls = line.lstyle,
                label = line.label,
                linewidth = line.width)

    # to avoid line overlapped on marker
    
    for i,line in enumerate( lines ):
        interpolator = create_linear_interpolator(line.df['x'], line.df['Cp'])
        ax.plot( x_sep[i],   interpolator(x_sep[i]),   'p', color=line.color, ms=15)
        ax.plot( x_att[i],   interpolator(x_att[i]),   'p', color=line.color, ms=15)
        ax.plot( x_pfmax[i], interpolator(x_pfmax[i]), '*', color=line.color, ms=20)

    ax.set_ylim([0.8,2.5])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.set_ylabel(r"$\langle p_w \rangle / p_{\infty}$")
    ax.tick_params(axis='y', pad=10)

    if show_label:
        ax.legend( fontsize=30 ) 
        
    adjust_plotting( ax )
    
    figname = "wall_pressure_awall"
    plt.savefig( figname + fmt )
    print( f"{figname.ljust(30)} is output in {OutPath}." )
    plt.show()


# ----------------------------------------------------------------------
# >>> plot pressure gradient                                     (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/01/06  - created
#
# Desc
#
# ----------------------------------------------------------------------


if plt_pwg:

    fig, ax = plt.subplots( figsize=figsize, constrained_layout=True )

    for i,line in enumerate( lines ):
        
        pwg = np.gradient( line.df['Cp'], line.df['x'] )
        
        ax.plot( line.df['x'], 
                 pwg,
                 line.color,
                 ls = line.lstyle,
                 label = line.label,
                 linewidth = line.width)

    # to avoid line overlapped on marker
    
    for i,line in enumerate( lines ):
        
        pwg = np.gradient( line.df['Cp'], line.df['x'] )
        interpolator = create_linear_interpolator(line.df['x'], pwg)
        ax.plot( x_sep[i],   interpolator(x_sep[i]),   'p', color=line.color, ms=15)
        ax.plot( x_att[i],   interpolator(x_att[i]),   'p', color=line.color, ms=15)
        ax.plot( x_pfmax[i], interpolator(x_pfmax[i]), '*', color=line.color, ms=20)

    ax.set_ylim([-0.1,0.6])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.set_ylabel(r"$\frac{d\langle p_w \rangle}{dx}/\frac{p_{\infty}}{\delta_0}$")
    ax.tick_params(axis='y', pad=10)

    if show_label:
        ax.legend( fontsize=30 ) 
        
    adjust_plotting( ax )

    figname = "wall_pressure_gradient"
        
    plt.savefig( figname + fmt )
    print( f"{figname.ljust(30)} is output in {OutPath}." )
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
    
    fig, ax = plt.subplots( figsize=figsize, constrained_layout=True )

    for i,line in enumerate( lines ):
        
        ax.plot( line.df['x'], 
                line.df['Cf'],
                line.color,
                ls = line.lstyle,
                label = line.label,
                linewidth = line.width)

    # to avoid line overlapped on marker
    
    for i,line in enumerate( lines ):
        interpolator = create_linear_interpolator(line.df['x'], line.df['Cf'])
        ax.plot( x_pfmax[i], interpolator(x_pfmax[i]), '*', color=line.color, ms=20)
            
    ax.plot( [-20,12],
             [0,0],
             'black',
             ls = '--',
             linewidth=2 )

    ax.set_ylim([-2.5,4.5])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    ax.set_ylabel("$C_fx10^3$")
    ax.tick_params(axis='y', pad=10)

    if show_label:
        ax.legend( fontsize=30 ) 
    
    adjust_plotting( ax )
        
    figname = "friction_coefficient_awall"
    
    plt.savefig( figname + fmt )
    print( f"{figname.ljust(30)} is output in {OutPath}." )
    plt.show()