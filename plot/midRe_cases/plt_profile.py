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
from   vista.tools        import find_indices

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 40

# =============================================================================

plt_u_vd           = False
plt_u              = False
plt_RS_uu          = False
plt_RS_vv          = False
plt_RS_ww          = False
plt_RS_uv          = False
plt_combined_RS    = False
plt_combined_RS_sn = False
plt_rho            = False
plt_T              = False
plt_Mt             = False
plt_Tt             = False
plt_pt             = False
plt_tke            = False
plt_p_fluc         = False
compute_p_fluc_w   = True    # compute the normalized pressure fluctuation

pure = False

location = 'upstream'   # 'upstream' or 'incip' 

fmt = '.png'

# =============================================================================

OutPath  = f'/home/wencan/temp/DataPost/midRe/profile_{location}/'

data0 = f'/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/profile_{location}'
data1 = f'/media/wencan/Expansion/temp/smooth_mid/postprocess/statistics/profile_{location}'
data2 = f'/media/wencan/Expansion/temp/220927/postprocess/statistics/profile_{location}'
data3 = f'/media/wencan/Expansion/temp/241030/postprocess/statistics/profile_{location}'
data4 = f'/media/wencan/Expansion/temp/231124/postprocess/statistics/profile_{location}'
data5 = f'/media/wencan/Expansion/temp/241018/postprocess/statistics/profile_{location}'

#data0 = '/home/wencanwu/my_simulation/temp/smooth_wall/x_-53.6.dat'

dataDNS = source_dir + '/database/Pirozzoli/M2_Retau_250'

datalist = [data0,             data2,             data1,             data4             ,data3                  ]
dy       = [0.0,               0.312,             0.0,               0.312             ,0.07921                ]
color    = ['gray',            'orangered',       'black',           'steelblue'       ,'yellowgreen'          ]
label    = [r'$\mathcal{LS}$', r'$\mathcal{LR}$', r'$\mathcal{HS}$', r'$\mathcal{HR}1$',r'$\mathcal{HR}2$'     ]
lstyle   = ['-',               '-.',              '--',              ':'               ,(0, (3, 1, 1, 1, 1, 1))]
width    = [4.0,               4.0,               4.0,               4.0               ,4.0                    ]
deltas   = [6.6732,            6.2661,            6.5022,            6.2369,           6.4216]
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
    
    line.label  = label[i]
    line.color  = color[i]
    line.width  = width[i]
    line.lstyle = lstyle[i]
    line.delta  = deltas[i]
    
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
        ax.legend( fontsize=30,frameon=False ) 

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
    y_lim         = [[-1,10], [-0.1,1.5], [-0.1,2.4], [-1,0.2]]
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
    y_lim         = [[-1,10], [-0.1,1.8], [-0.1,2.5], [-1.4,0.2]]
    major_locator = [2.0, 0.4, 0.5, 0.4]
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
    

# ----------------------------------------------------------------------
# >>> plot total temperature                                     (Nr.)
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

if plt_Tt:
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['Tt'],
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
    
    figname = 'Tt'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.set_ylabel( r'$T_{t}$' )
        ax.tick_params( axis='x', pad = 15 )
        ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()

# ----------------------------------------------------------------------
# >>> plot total pressure                                        (Nr.)
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

if plt_pt:
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['pt'],
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
    
    figname = 'pt'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
        ax.set_ylabel( r'$p_{t}$' )
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

# ----------------------------------------------------------------------
# >>> plot p`                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/05/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_p_fluc:

    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )

    for line in lines:
        
        p_fluc = np.array(line.df['p`'])
        
        ax.plot( line.df['ys+'], 
                 (p_fluc/line.tau_ave)**2,
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,20000] )
    ax.set_ylim( [0.0,15]   )

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
    
    figname = 'p_fluc_squared'
    
    if pure:    
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$",labelpad=-5 )  
        ax.set_ylabel( r"$\langle p'p' \rangle / {\tau_w}^2$" )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()
    
if compute_p_fluc_w:
    
    for line in lines:
        
        il, ir = find_indices(line.df['ys+'],15.0, mode='sequential')
        
        rho15 = line.df['rho'][il] + (line.df['rho'][ir]-line.df['rho'][il]) * (15.0-line.df['ys+'][il])/(line.df['ys+'][ir]-line.df['ys+'][il])
        mu15  = line.df['mu'][il] + (line.df['mu'][ir]-line.df['mu'][il]) * (15.0-line.df['ys+'][il])/(line.df['ys+'][ir]-line.df['ys+'][il])
        
        u_taw15 = np.sqrt(line.tau_ave/rho15)
        
        Re15 = rho15 * u_taw15 * line.delta / mu15
        
        T_w = line.df['T'][0]
        c_w = np.sqrt(1.4 * 287.0508571 * T_w)
        
        mach_tau = line.u_tau / c_w
        
        pp_w_plus = 20.25 + \
                    (-87.3*Re15**(-0.25) + 94.09*Re15**(-0.5)) + \
                    2.4 * mach_tau ** 2  + 8312.5 * mach_tau ** 4
        
        print(f"line {line.label} : p'p' = {pp_w_plus:.2f}")
        
        
        