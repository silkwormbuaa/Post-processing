#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   to_reviewer_profile.py
@Time    :   2025/10/15 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plot the profiles at some local points, like ridge/valley at 
             the upstream and incipient interaction location.
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
from   vista.directories  import create_folder

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 40

# =============================================================================

plt_u_vd           = True
plt_u              = True
plt_RS_uu          = True
plt_RS_vv          = True
plt_RS_ww          = True
plt_RS_uv          = True
plt_combined_RS    = True
plt_combined_RS_sn = True
plt_T              = True
plt_tke            = True

pure  = False

fmt = '.png'

# =============================================================================

OutPath  = '/home/wencan/temp/DataPost/midRe/profile_local/upstream_plane/jfm_reviewer'

cases = ['smooth_adiabatic', '220927', '220927','smooth_mid', '231124','231124','241030','241030']

#dy       = [0.0,               0.312,             0.0,               0.312             ,0.07921                ]
dy       = [0.0,               0.0,              0.52,               0.0,               0.0               ,0.52              ,0.0               ,0.13679           ]
color    = ['gray',            'orangered',       'orangered',       'black',           'steelblue'       ,'steelblue'       ,'yellowgreen'     ,'yellowgreen'     ]
label    = [r'$\mathcal{LS}$', r'$\mathcal{LR}_r$', r'$\mathcal{LR}_v$', r'$\mathcal{HS}$', r'$\mathcal{HR}1_r$',r'$\mathcal{HR}1_v$',r'$\mathcal{HR}2_r$',r'$\mathcal{HR}2_v$']
locs     = ['ridge',          'ridge',           'valley',           'ridge',           'ridge'           ,'valley'          ,     'ridge'      ,'valley'          ]
lstyle   = ['-',               '-',               '--',              '--',              '-'               ,'--'              , '-'              , '--'             ]
width    = [4.0,               4.0,               4.0,               4.0,               4.0               ,4.0               ,4.0               ,4.0               ]
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

for i, datapath in enumerate(cases):
    
    os.chdir('/home/wencan/temp/'+datapath+'/postprocess/statistics/profile_upstream')
    
    if locs[i] == 'ridge':
        datafile = 'profile_upstream_ridge.dat'
        line = ProfileData( datafile )
    else:
        datafile = 'profile_upstream_valley.dat'
        line = ProfileData( datafile )
    
    line.shift_y( dy[i] )
    line.inner_normalize('wall_statistics.dat')
    line.vd_transform()
    
    line.label = label[i]
    line.color = color[i]
    line.width = width[i]
    line.lstyle = lstyle[i]
    
    if locs[i] == 'ridge':
        normalized_file = 'profile_normalized_ridge.dat'
    else:
        normalized_file = 'profile_normalized_valley.dat'
    
    line.df.to_string( normalized_file,
                       index=False, 
                       float_format='%15.7f',
                       justify='left' )
    
    lines.append(line)


os.chdir( create_folder(OutPath) )


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
    
    figname = 'u_vd_legend'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.tick_params( axis='x', pad=15 )
        ax.set_ylabel( r'$\langle u \rangle ^+_{vD}$' )
        ax.tick_params( axis='y', pad=10 )
        ax.legend( fontsize=20,frameon=False, loc='lower right' ) 

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
    ax.set_ylim( [-1,12]   )

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
    ax.set_ylim( [-0.1,1.8]   )

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
    ax.set_ylim( [-0.1,2.7]   )

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
    ax.set_ylim( [-1.5,0.2] )

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
# >>> plot combined Reynolds stresses                          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/03/05  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_combined_RS:
    
    fig = plt.figure(figsize=[20,12], constrained_layout=True)
    gs  = fig.add_gridspec(2, 2)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    
    axs           = [ax1, ax2, ax3, ax4]
    vars          = ['u`u`+', 'v`v`+', 'w`w`+', 'u`v`+']
    y_lim         = [[-1,12], [-0.1,1.8], [-0.1,2.7], [-1.5,0.2]]
    major_locator = [2.0, 0.4, 0.5, 0.4]
    ylabels       = [r"$\langle \rho \rangle \langle u^{'} u^{'} \rangle / \tau_w$",
                     r"$\langle \rho \rangle \langle v^{'} v^{'} \rangle / \tau_w$",
                     r"$\langle \rho \rangle \langle w^{'} w^{'} \rangle / \tau_w$",
                     r"$\langle \rho \rangle \langle u^{'} v^{'} \rangle / \tau_w$"]
    
    for i, ax in enumerate( axs ):
    
        for line in lines:
            ax.plot( line.df['ys+'],
                     line.df[vars[i]],
                     line.color,
                     label = line.label,
                     ls    = line.lstyle,
                     linewidth = line.width)
    
        ax.set_xscale( "symlog", linthresh=1 )
        ax.set_xlim( [1,3000] )
        ax.set_ylim( y_lim[i] )
        ax.minorticks_on()
        ax.tick_params( which='major',
                        axis='both',
                        direction='in',
                        length=20,
                        width=1.5,
                        pad=10)
        ax.tick_params( which='minor',
                        axis='both', 
                        direction='in',
                        length=10,
                        width=1.5)
        if i == 0 or i == 1:
            ax.set_xticklabels([])
        x_minor = matplotlib.ticker.LogLocator( 
                            base = 10.0, subs = np.arange(1.0,10.0), numticks=100 )
        ax.xaxis.set_minor_locator( x_minor )
        ax.yaxis.set_major_locator(ticker.MultipleLocator(major_locator[i]))
        
        ax.set_ylabel( ylabels[i] )
        
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(3)
    
    ax1.legend( frameon=False, loc='upper right', fontsize=25 )
    
    ax3.set_xlabel( "$y_s^+$", labelpad=-5 )
    ax4.set_xlabel( "$y_s^+$", labelpad=-5 )
        
    plt.savefig( "RS_combined_legend" + fmt )
    plt.show()
    

# ----------------------------------------------------------------------
# >>> plot combined Reynolds stresses (normalize with smooth wall) (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/03/05  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_combined_RS_sn:
    
    fig = plt.figure(figsize=[20,12], constrained_layout=True)
    gs  = fig.add_gridspec(2, 2)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    
    axs           = [ax1, ax2, ax3, ax4]
    vars          = ['u`u`+', 'v`v`+', 'w`w`+', 'u`v`+']
    y_lim         = [[-1,15], [-0.1,2.2], [-0.1,3.5], [-2.2,0.2]]
    major_locator = [4.0, 0.5, 1.0, 0.5]
    ylabels       = [r"$\langle \rho \rangle \langle u^{'} u^{'} \rangle / \tau_{sw}$",
                     r"$\langle \rho \rangle \langle v^{'} v^{'} \rangle / \tau_{sw}$",
                     r"$\langle \rho \rangle \langle w^{'} w^{'} \rangle / \tau_{sw}$",
                     r"$\langle \rho \rangle \langle u^{'} v^{'} \rangle / \tau_{sw}$"]
    
    for i, ax in enumerate( axs ):
    
        for j, line in enumerate(lines):
            
            if j in [0,1]:
                value = line.df[vars[i]]* line.tau_ave / lines[0].tau_ave
            if j in [2,3,4]:
                value = line.df[vars[i]]* line.tau_ave / lines[2].tau_ave
            
            ax.plot( line.df['ys+'],
                     value,
                     line.color,
                     label = line.label,
                     ls    = line.lstyle,
                     linewidth = line.width)
    
        ax.set_xscale( "symlog", linthresh=1 )
        ax.set_xlim( [1,3000] )
        ax.set_ylim( y_lim[i] )
        ax.minorticks_on()
        ax.tick_params( which='major',
                        axis='both',
                        direction='in',
                        length=20,
                        width=1.5,
                        pad=10)
        ax.tick_params( which='minor',
                        axis='both', 
                        direction='in',
                        length=10,
                        width=1.5)
        if i == 0 or i == 1:
            ax.set_xticklabels([])
        x_minor = matplotlib.ticker.LogLocator( 
                            base = 10.0, subs = np.arange(1.0,10.0), numticks=100 )
        ax.xaxis.set_minor_locator( x_minor )
        ax.yaxis.set_major_locator(ticker.MultipleLocator(major_locator[i]))
        
        ax.set_ylabel( ylabels[i] )
        
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(3)
    
    ax1.legend( frameon=False, loc='upper right', fontsize=25 )
    
    ax3.set_xlabel( "$y_s^+$", labelpad=-5 )
    ax4.set_xlabel( "$y_s^+$", labelpad=-5 )
        
    plt.savefig( "RS_combined_sn" + fmt )
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
    
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)

    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['T'],
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
#    ax.set_ylim( [0,25] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0), numticks=100 )
    
    ax.xaxis.set_minor_locator( x_minor )


    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'T'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.set_ylabel( r'$T$' )
        ax.tick_params( axis='x', pad = 15 )
        ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> plot tke                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/04/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_tke:

    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )

    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u`u`+'] + line.df['v`v`+'] + line.df['w`w`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,3000] )
    ax.set_ylim( [-1,15]   )

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
    
    figname = 'tke'
    
    if pure:    
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$",labelpad=-5 )  
        ax.set_ylabel( r"$\rho \langle tke \rangle / \tau_w$" )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y',pad = 10)
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()

