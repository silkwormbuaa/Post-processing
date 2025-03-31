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
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import LineData
from   vista.tools       import create_linear_interpolator
from   vista.directories import create_folder

plt.rcParams["text.usetex"]         = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family']         = "Times New Roman"
plt.rcParams['font.size']           = 40

OutPath  = "/home/wencan/temp/DataPost/midRe/averaged_streamwise_vars/with_241028"

data0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data1 = '/media/wencan/Expansion/temp/smooth_mid/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data2 = '/media/wencan/Expansion/temp/220927/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data3 = '/media/wencan/Expansion/temp/241030/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data4 = '/media/wencan/Expansion/temp/231124/postprocess/statistics/wall_projection/streamwise_vars.pkl'
data5 = '/media/wencan/Expansion/temp/241018/postprocess/statistics/wall_projection/streamwise_vars.pkl'

sw_iw_pw_file = '/media/wencan/Expansion/temp/data_luis/midRe_streamwise/pw_SBLI_B2.dat'
sw_iw_pwrms_file = '/media/wencan/Expansion/temp/data_luis/midRe_streamwise/pw_rms_SBLI_B2.dat'
sw_iw_Cf_file = '/media/wencan/Expansion/temp/data_luis/midRe_streamwise/Cf_SBLI_B2.dat'

add_sw_luis = False

# color, label, lstyle, width's order will be adjusted automatically
# x_sep, x_att, x_pfmax should be carefully checked. 

datalist = [data0,              data2,             data1,             data4,              data3                  ,data5]
color    = ['gray',             'orangered',       'black',           'steelblue',        'yellowgreen'          ,'red'] 
label    = [r'$\mathcal{LS}$',  r'$\mathcal{LR}$', r'$\mathcal{HS}$', r'$\mathcal{HR}1$', r'$\mathcal{HR}2$'     ,r'$\mathcal{HR}3$']
lstyle   = ['-',                '-.',               '--',             ':',                (0, (3, 1, 1, 1, 1, 1)),'--']
width    = [4.0,                 4.0,               4.0,                 4.0,                4.0                 ,4.0]
x_sep    = [-8.42,              -10.84,            -7.33,              -9.96,             -10.36                 ,-15.03]
x_att    = [1.10,                 2.45,             1.89,               2.86,               2.56                 ,3.57]
x_pfmax  = [-7.136,             -10.06,            -8.53,             -11.58,             -10.16                 ,-14.75]

lines    = []

plt_pwfluc     = True     # pressure fluctuation
plt_pwfluc_ln  = True     # locally normalized
plt_pw         = True     # wall pressure
plt_pwg        = True     # pressure gradient
plt_Cf         = True

pure       = False
show_label = False

figsize    = [15,8]
fmt        =  '.png' # or '.png'

# - read in data files

for i, datafile in enumerate( datalist ):
    
    line = LineData()
    
    with open(datafile,'rb') as f:  line.df = pickle.load( f )
    line.color  = color[i]
    line.label  = label[i]
    line.lstyle = lstyle[i]
    line.width  = width[i]
    
    lines.append( line )

os.chdir( create_folder(OutPath) )

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
# 2024/08/16  - created
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
        ax.plot( x_sep[i],   interpolator(x_sep[i]),   'p', color=line.color, ms=15)
        ax.plot( x_att[i],   interpolator(x_att[i]),   'p', color=line.color, ms=15)
        ax.plot( x_pfmax[i], interpolator(x_pfmax[i]), '*', color=line.color, ms=15)

    if add_sw_luis:
        swdf = pd.read_csv(sw_iw_pwrms_file, delimiter=r'\s+')
        ax.plot( swdf['x']*7.15/5.2, 
                 swdf['y'],
                 'black',
                 ls = '-',
                 label = 'midRe_smooth',
                 linewidth = 4.0)


    figname = "pressure_fluctuation"

    # adjust the y axis

    ax.set_ylim([0.01,0.10])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
    ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_{\infty}$" )
    ax.tick_params(axis='y', pad=10)
    
    adjust_plotting( ax )
    
    if show_label:
        ax.legend( fontsize=20 ) 
        
    plt.savefig( figname + fmt )
    print( f"{figname.ljust(25)} is output in {OutPath}." )
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

    if add_sw_luis:
        swdf = pd.read_csv(sw_iw_pwrms_file, delimiter=r'\s+')
        ax.plot( swdf['x']*7.15/5.2, 
                 swdf['y'],
                 'black',
                 ls = '-',
                 label = 'midRe_smooth',
                 linewidth = 4.0)

    figname = "pressure_fluctuation_ln"

    # adjust the y axis

    ax.set_ylim([0.01,0.08])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
    ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_{w}$" )
    ax.tick_params(axis='y', pad=10)
    
    adjust_plotting( ax )
    
    if show_label:
        ax.legend( fontsize=20 ) 
        
    plt.savefig( figname + fmt )
    print( f"{figname.ljust(25)} is output in {OutPath}." )
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
        
    if add_sw_luis:
        swdf = pd.read_csv(sw_iw_pw_file, delimiter=r'\s+')
        ax.plot( swdf['x']*7.15/5.2, 
                 swdf['y'],
                 'black',
                 ls = '-',
                 label = 'midRe_smooth',
                 linewidth = 4.0)
        
    ax.set_ylim([0.8,2.5])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    
    ax.set_ylabel(r"$\langle p_w \rangle/p_{\infty}$")
    ax.tick_params(axis='y', pad=10)

    if show_label:
        ax.legend( fontsize=30 ) 

    adjust_plotting( ax )

    figname = "wall_pressure"

    plt.savefig( figname + fmt )
    print( f"{figname.ljust(25)} is output in {OutPath}." )
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

    ax.set_ylabel(r"$\frac{d \langle p_w \rangle }{dx}/p_{\infty}$")
    ax.tick_params(axis='y', pad=10)

    if show_label:
        ax.legend( fontsize=30 ) 

    adjust_plotting( ax )

    figname = "wall_pressure_gradient"
    
    plt.savefig( figname + fmt )
    print( f"{figname.ljust(25)} is output in {OutPath}." )
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
    
    ax.set_ylim([-2.5,5.0])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    
    ax.set_ylabel(r"$C_f \times 10^3$")
    ax.tick_params(axis='y', pad=10)

    if show_label:
        ax.legend( fontsize=30 ) 

    adjust_plotting( ax )

    figname = "friction_coefficient"
        
    plt.savefig( figname + fmt )
    print( f"{figname.ljust(25)} is output in {OutPath}." )
    plt.show()