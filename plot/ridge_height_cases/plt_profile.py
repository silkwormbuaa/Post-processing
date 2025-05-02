#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_profile.py
@Time    :   2024/03/11 
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

plt_u_vd_lw = True
plt_u_vd    = True
plt_u       = True
plt_RS_uu   = True
plt_RS_vv   = True
plt_RS_ww   = True
plt_RS_uv   = True
plt_RS_DNS  = True
plt_rho     = True
plt_T       = False
plt_Mt      = False

pure = False

fmt = '.png'

# =============================================================================

OutPath  = '/media/wencan/Expansion/temp/DataPost/lowRe_ridge_height/profiles_upstream/DLES14_black_white'

data0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/profile_upstream'
data1 = '/media/wencan/Expansion/temp/240211/postprocess/statistics/profile_upstream'
data2 = '/media/wencan/Expansion/temp/220927/postprocess/statistics/profile_upstream'
data3 = '/media/wencan/Expansion/temp/240210/postprocess/statistics/profile_upstream'


#data0 = '/home/wencanwu/my_simulation/temp/smooth_wall/x_-53.6.dat'

dataDNS = source_dir + '/database/Pirozzoli/M2_Retau_250'

datalist = [data0, data1,   data2,   data3]
dy       = [0,     0.156,   0.312,   0.624]
color    = ['black','black', 'black',  'black']
label    = ['',    '0.05',   '0.1',  '0.20']
lstyle   = ['--',  ':',     '-.',    (0, (3, 1, 1, 1, 1, 1)) ]
width    = [4.0,   4.0,      4.0,    4.0 ]
lines = []

# ----------------------------------------------------------------------
# >>> Initialize data                                            ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/11  - created
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
    
    line.label = r'$H/\delta_0=$' + label[i]
    line.color = color[i]
    line.width = width[i]
    line.lstyle = lstyle[i]
    
    line.df.to_string( "profile_normalized.dat",
                       index=False, 
                       float_format='%15.7f',
                       justify='left' )
    
    lines.append(line)
    
lines[0].label = 'smooth'


lineDNS = ProfileData( dataDNS )
lineDNS.df = lineDNS.sparse_log('y+', 0.04)
lineDNS.label = "DNS"
lineDNS.color = 'black'
lineDNS.lstyle = 's'
lineDNS.marker = markers.MarkerStyle(marker='s')

os.chdir( OutPath )

# ----------------------------------------------------------------------
# >>> plot comparison between smooth wall and law of wall        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_u_vd_lw:
    
    fig, ax = plt.subplots(figsize=[9,8], constrained_layout=True)

    x1 = np.linspace(1,12,50)
    y1 = np.linspace(1,12,50)
    
    x2 = np.linspace(8,400,300)
    y2  = np.log(x2)/0.41 + 5.1    
    
    ax.plot( x1,y1,'black',ls=':')
    ax.plot( x2,y2,'black',ls=':')
       
    ax.plot( lineDNS.df['y+'], 
             lineDNS.df['u_vd+'],
             lineDNS.color, 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')

    ax.plot( lines[0].df['y+'], 
             lines[0].df['u+_vd'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)
    
    ax.text( 1.8, 8, r'$u^+=y^+$', fontsize=35 )
    ax.text( 1.7, 17,r"$u^+=\frac{1}{0.41}ln(y^+)+5.1$",fontsize=35 )
    
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
    
    ax.set_xlim( [1,1000] )
    ax.set_ylim( [0,23] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )

#    ax.grid(visible=False, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)
    
    figname = 'Uvd_DNS'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel(  r"$y_s^+$",labelpad=-5 )  
        ax.tick_params( axis='x', pad=15 )
        ax.set_ylabel( r'$\langle u \rangle ^+_{vd}$' )
        ax.tick_params( axis='y', pad=10 )
#        ax.legend( prop={'size':30} ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)

    plt.savefig( figname + fmt )
    plt.show()    



# ----------------------------------------------------------------------
# >>> Plot van Driest transformed u profile                      ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/11  - created
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
    
    ax.set_xlim( [1,1000] )
    ax.set_ylim( [0,23] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )


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
#        ax.legend( ) 

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
# 2024/03/11  - created
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
    
    ax.set_xlim( [1,1000] )
    ax.set_ylim( [0,23] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
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
# 2024/03/11  - created
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


    ax.set_xlim( [1,1000] )
    ax.set_ylim( [-1,8]   )

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
                        base = 10.0, subs = np.arange(1.0,10.0) )
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


    ax.set_xlim( [1,1000] )
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
                        base = 10.0, subs = np.arange(1.0,10.0) )
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


    ax.set_xlim( [1,1000] )
    ax.set_ylim( [-0.1,2]   )

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
                        base = 10.0, subs = np.arange(1.0,10.0) )
    ax.xaxis.set_minor_locator( x_minor )
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))


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
# 2024/03/11  - created
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


    ax.set_xlim( [1,1000] )
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
                        base = 10.0, subs = np.arange(1.0,10.0) )
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
# >>> plot RS comparison between DNS and smooth wall             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_RS_DNS:
    
    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )

    ax.plot( lineDNS.df['y+'], 
             lineDNS.df['urms+']**2*lineDNS.df['sqrt(rho/rho_w)']**2,
             'black',
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')

    ax.plot( lineDNS.df['y+'], 
             lineDNS.df['vrms+']**2*lineDNS.df['sqrt(rho/rho_w)']**2,
             'black', 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')
    
    ax.plot( lineDNS.df['y+'], 
             (lineDNS.df['wrms+'])**2*lineDNS.df['sqrt(rho/rho_w)']**2,
             'black', 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')
    
    ax.plot( lineDNS.df['y+'], 
             (lineDNS.df['uv+'])*lineDNS.df['sqrt(rho/rho_w)']**2,
             'black', 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')
    
    ax.plot( lines[0].df['y+'], 
             lines[0].df['u`u`+'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)

    ax.plot( lines[0].df['y+'], 
             lines[0].df['v`v`+'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)

    ax.plot( lines[0].df['y+'], 
             lines[0].df['w`w`+'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)

    ax.plot( lines[0].df['y+'], 
             lines[0].df['u`v`+'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)
        
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
    
    ax.set_xlim( [1,1000] )
    ax.set_ylim( [-1.5,8.5] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )


#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_DNS'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel(  r"$y_s^+$",labelpad=-5 )  
        ax.tick_params( axis='x', pad=15 )
        ax.set_ylabel( r"$\rho \langle u_i^{'}u_j^{'}\rangle /\tau_w$" )
        ax.tick_params( axis='y', pad=10 )
#        ax.legend( prop={'size':22} ) 

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
# 2024/03/11  - created
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
    
    ax.set_xlim( [1,1000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
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
# 2024/03/11  - created
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
    
    ax.set_xlim( [1,1000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
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
# 2024/03/11  - created
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
    
    ax.set_xlim( [1,1000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )

#    ax.legend( ) 
#    ax.set_title( r"$M_t$ profile" )

    ax.grid()
    
    plt.savefig( "Mxt_profile_shifted" + fmt )
    plt.show()