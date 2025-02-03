#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_profile.py
@Time    :   2024/08/16 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plot the profiles of the low Re different ridge cases
'''

import os
import sys

import matplotlib
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker
import matplotlib.markers as     markers
import numpy              as     np

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line         import ProfileData

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 40

# =============================================================================

plt_u_vd    = True
plt_u       = True
plt_RS_uu   = True
plt_RS_vv   = True
plt_RS_ww   = True
plt_RS_uv   = True
plt_rho     = True
plt_T       = False
plt_Mt      = False

pure = False

fmt = '.png'

# =============================================================================

OutPath  = '/home/wencan/temp/DataPost/midRe/upstream_profile/'

data0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/upstream_profile'
data1 = '/media/wencan/Expansion/temp/smooth_mid/postprocess/statistics/upstream_profile'
data2 = '/media/wencan/Expansion/temp/220927/postprocess/statistics/upstream_profile'
data3 = '/media/wencan/Expansion/temp/241030/postprocess/statistics/upstream_profile'
data4 = '/media/wencan/Expansion/temp/231124/postprocess/statistics/upstream_profile'
data5 = '/media/wencan/Expansion/temp/241018/postprocess/statistics/upstream_profile'

#data0 = '/home/wencanwu/my_simulation/temp/smooth_wall/x_-53.6.dat'

dataDNS = source_dir + '/database/Pirozzoli/M2_Retau_250'

datalist = [data0,          data1,         data2,   data3,   data4]
dy       = [0.0,            0.0,           0.312,   0.07921, 0.312]
color    = ['gray',         'black',       'orangered',  'yellowgreen', 'steelblue']
label    = ['lowRe_smooth', 'midRe_smoth', 'lowRe_rough','midRe_0.026', 'midRe_0.1']
lstyle   = ['-',            '--',             '-.',      (0, (3, 1, 1, 1, 1, 1)),  ':']
width    = [4.0,            4.0,           4.0,        4.0,        4.0]
lines = []

# ----------------------------------------------------------------------
# >>> Initialize data                                            ( 1 )
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

for i, datapath in enumerate(datalist):
    
    os.chdir(datapath)
    
    line = ProfileData('profile_mean.dat') 
    line.shift_y( dy[i] )
    line.inner_normalize('wall_statistics.dat')
    line.vd_transform()
    
    line.label = label[i]
    line.color = color[i]
    line.width = width[i]
    line.lstyle = lstyle[i]
    
    line.df.to_string( "profile_normalized.dat",
                       index=False, 
                       float_format='%15.7f',
                       justify='left' )
    
    lines.append(line)
    

lineDNS = ProfileData( dataDNS )
lineDNS.df = lineDNS.sparse_log('y+', 0.04)
lineDNS.label = "DNS"
lineDNS.color = 'black'
lineDNS.lstyle = 's'
lineDNS.marker = markers.MarkerStyle(marker='s')

print( OutPath)
os.chdir( OutPath )


# ----------------------------------------------------------------------
# >>> Plot van Driest transformed u profile                      ( 1 )
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

if plt_u_vd :
    
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u+_vd'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )

    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=1.5)
    
    ax.set_xlim( [1,3000] )
    ax.set_ylim( [0,26] )
    
    x_minor = ticker.LogLocator( base=10.0, subs = np.arange(1.0,10.0), numticks=100 )
    
    ax.xaxis.set_minor_locator( x_minor )
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'Uvd'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.tick_params( axis='x', pad=15 )
        ax.set_ylabel( r'$\langle u \rangle ^+_{vd}$' )
        ax.tick_params( axis='y', pad=10 )
#        ax.legend( fontsize=30 ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)

    plt.savefig( figname + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot u profile                                             ( 2 )
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

if plt_u :
    
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )

    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.5)
    
    ax.set_xlim( [1,3000] )
    ax.set_ylim( [0,25] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0), numticks=100 )
    
    ax.xaxis.set_minor_locator( x_minor )


    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'u'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.set_ylabel( r'$u^+$' )
        ax.tick_params( axis='x', pad = 15 )
        ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot Reynolds Stress Profile                              ( 3 )
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

if plt_RS_uu:

    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )

    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u`u`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,3000] )
    ax.set_ylim( [-1,10]   )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0), numticks=100 )
    ax.xaxis.set_minor_locator( x_minor )

    # set spacing between major tickers.
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))

#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_uu'
    
    if pure:    
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$",labelpad=-5 )  
        ax.set_ylabel( r"$\rho \langle u^{'} u^{'} \rangle / \tau_w$" )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y',pad = 10)
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


if plt_RS_vv:
    
    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['v`v`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
    
    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,3000] )
    ax.set_ylim( [-0.1,1.5]   )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0), numticks=100 )
    ax.xaxis.set_minor_locator( x_minor )
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))


#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_vv'
    
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.set_ylabel( r"$\rho \langle v^{'} v^{'} \rangle / \tau_w$" )
        ax.tick_params( axis='x', pad = 15 )
        ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()
    

if plt_RS_ww:
    
    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['w`w`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
    
    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,3000] )
    ax.set_ylim( [-0.1,2.4]   )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0), numticks=100 )
    ax.xaxis.set_minor_locator( x_minor )
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))


#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_ww'
    
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.set_ylabel( r"$\rho \langle w^{'} w^{'} \rangle / \tau_w$")
        ax.tick_params( axis='x', pad = 15 )
        ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> plot Reynolds stress                                        (Nr.)
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

if plt_RS_uv:

    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u`v`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
    
    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,3000] )
    ax.set_ylim( [-1,0.2] )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0), numticks=100 )
    ax.xaxis.set_minor_locator( x_minor )
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))

#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_uv'
    
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.set_ylabel( r"$\rho \langle u^{'} v^{'} \rangle / \tau_w$" )
        ax.tick_params( axis='x', pad = 15 )
        ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()
    

# ----------------------------------------------------------------------
# >>> Plot rho profile                                           ( 4 )
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

if plt_rho :
    
    fig, ax = plt.subplots(figsize=[8,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['rho'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$" )  
    ax.tick_params( axis='x' )
    
    ax.set_ylabel( r"$\rho$" )
    ax.tick_params( axis='y' )
    
    ax.set_xlim( [1,3000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0), numticks=100 )
    
    ax.xaxis.set_minor_locator( x_minor )

#    ax.legend( ) 
    ax.set_title( r"$\rho$ profile" )
    
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    ax.grid()
    
    plt.savefig( "rho_shifted" + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot temperature profile                                  ( 5 )
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

if plt_T :
    
    fig, ax = plt.subplots(figsize=[8,8])

    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['T'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
        
    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$" )  
    ax.tick_params( axis='x' )
    
    ax.set_ylabel( r'$T/T_{\infty}$' )
    ax.tick_params( axis='y' )
    
    ax.set_xlim( [1,3000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0), numticks=100 )
    
    ax.xaxis.set_minor_locator( x_minor )

#    ax.legend( ) 
#    ax.set_title( r"$u^+_{VD}$ profile", size=20 )

    ax.grid()
    
    plt.savefig( "T_shifted_new" + fmt)
    plt.show()

# ----------------------------------------------------------------------
# >>> Plot Mach_prime(fluctuating Mach) number                   ( 6 )
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

if plt_Mt :
    
    fig, ax = plt.subplots(figsize=[8,8])

    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$" )  
    ax.tick_params( axis='x' )
    
    ax.set_ylabel( r'$M_t$' )
    ax.tick_params( axis='y' )
    
    ax.set_xlim( [1,3000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0), numticks=100 )
    
    ax.xaxis.set_minor_locator( x_minor )

#    ax.legend( ) 
#    ax.set_title( r"$M_t$ profile" )

    ax.grid()
    
    plt.savefig( "Mxt_profile_shifted" + fmt )
    plt.show()