# -*- coding: utf-8 -*-
'''
@File    :   plt_streamwise_vars.py
@Time    :   2025/07/01 
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
from   vista.params      import Params
from   vista.directories import Directories
from   vista.tools       import create_linear_interpolator
from   vista.directories import create_folder


plt.rcParams["text.usetex"]         = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family']         = "Times New Roman"
plt.rcParams['font.size']           = 40

OutPath  = "/media/wencan/Expansion/temp/DataPost/herringbones_patch/averaged_streamwise_vars_h10_h03"

# sw_pfluc_file = '/home/wencanwu/my_simulation/temp/smooth_wall/line_p_prime.dat'
# sw_Cf_file = '/home/wencanwu/my_simulation/temp/smooth_wall/x_cf_STBLI_Wencan.dat'
# sw_Cp_file = '/home/wencanwu/my_simulation/temp/smooth_wall/Cf_flat_new.dat'

cases = ['smooth_adiabatic','250710','250821']

casepaths = [f'/media/wencan/Expansion/temp/{case}' for case in cases]
datapaths = [f'{casepath}/postprocess/statistics/wall_projection/streamwise_vars.pkl' for casepath in casepaths]

dy       = [0.0,       0.0, 0.0]
color    = ['gray', 'blue', 'red']
label    = ['smooth',   'cd1', 'cd2']
lstyle   = [':',     '-.', '--']
width    = [4.0,      4.0, 4.0]

x_sep, x_att, x_pfmax = [], [], []

for casepath in casepaths:
    dirs   = Directories( casepath )
    params = Params( dirs.case_para_file )
    x_sep.append(params.x_sep)
    x_att.append(params.x_att)
    x_pfmax.append(params.x_pfmax) 

lines    = []

plt_pwfluc    = True
plt_pwfluc_ln = True
plt_pw        = True
plt_pwg       = True
plt_Cf        = True

pure          = False
show_label    = False

figsize       = [15,8]
fmt           =  '.png' # or '.png'

# - read in data files

for i, datafile in enumerate( datapaths ):
    
    line = LineData()
    
    with open(datafile,'rb') as f:  line.df = pickle.load( f )
    line.color  = color[i]
    line.label  = label[i]
    line.lstyle = lstyle[i]
    line.width  = width[i]
    
    lines.append( line )


os.chdir(create_folder(OutPath))

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