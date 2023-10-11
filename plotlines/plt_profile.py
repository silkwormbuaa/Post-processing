# -*- coding: utf-8 -*-
'''
@File    :   plt_profile.py
@Time    :   2022/10/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   T, and Mt profile not yet finished.
'''

import os
import sys

import matplotlib
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker
import matplotlib.markers as     markers

import numpy              as     np
import pandas             as     pd

source_dir = os.path.realpath(__file__).split('plotlines')[0]
sys.path.append( source_dir )

from   vista.line         import ProfileData


OutPath  = '/home/wencanwu/my_simulation/temp/DataPost/profiles'

data1 = '/media/wencanwu/Seagate Expansion Drive/temp/221014/results/profile'
data2 = '/media/wencanwu/Seagate Expansion Drive/temp/220926/results/profile'
data3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/results/profile'
data4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/profile'
data5 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results/profile'

data0 = '/home/wencanwu/my_simulation/temp/smooth_wall/x_-53.6.dat'

dataDNS = source_dir + '/database/Pirozzoli/M2_Retau_250'

datalist = [data1,   data2,   data3,                   data4,        data5]
dy       = [0.494,   0.468,   0.416,                   0.312,        0.26]
color    = ['green', 'blue', 'black',                  'red',        'purple']
label    = ['2.0',   '1.0',  '0.5',                    '0.25',       '0.125']
lstyle   = [':',     '-.',    (0, (3, 1, 1, 1, 1, 1)), (0, (10, 3)), '-']
width    = [4.0,      4.0,    4.0,                     4.0,          4.0]
lines = []

plt_u_vd_lw = False
plt_u_vd  = False
plt_u     = False
plt_RS_uu = False
plt_RS_vv = False
plt_RS_ww = False
plt_RS_uv = False
plt_rho   = False
plt_T     = False
plt_Mt    = False

pure = True

# ----------------------------------------------------------------------
# >>> Initialize data                                            ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/01  - created
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
    
    line.label = r'$\mathrm{D/\delta_0=}$' + label[i]
    line.color = color[i]
    line.width = width[i]
    line.lstyle = lstyle[i]
    
    line.df.to_string( "profile_normalized.dat",
                       index=False, 
                       float_format='%15.7f',
                       justify='left' )
    
    lines.append(line)
    

line0 = ProfileData(data0)
line0.rho_ave = line0.df['rho'][0]
line0.vd_transform()
line0.label = 'smooth'
line0.color = 'gray'
line0.width = 4.0
line0.lstyle = '--'

lineDNS = ProfileData( dataDNS )
lineDNS.label = "DNS"
lineDNS.color = 'gray'
lineDNS.lstyle = 's'
lineDNS.marker = markers.MarkerStyle(marker='s')

os.chdir( OutPath )

# ----------------------------------------------------------------------
# >>> plot comparison between smooth wall and law of wall                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/04  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_u_vd_lw:
    
    x1 = np.linspace(1,12,50)
    y1 = np.linspace(1,12,50)
    
    x2 = np.linspace(8,400,300)
    y2  = np.log(x2)/0.41 + 5.1    

    fig, ax = plt.subplots(figsize=[10,8])
    
    ax.plot( x1,y1,'black',ls=':')
    ax.plot( x2,y2,'black',ls=':')
       
    ax.plot( line0.df['y+'], 
             line0.df['u+_vd'],
             line0.color, 
             label = line0.label, 
             ls    = line0.lstyle,
             linewidth = line0.width)

    ax.plot( lineDNS.df['y+'][::2], 
             lineDNS.df['u_vd+'][::2],
             lineDNS.color, 
             marker='s',
             fillstyle='none')
    
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
    
    figname = 'law_wall'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
        ax.tick_params( axis='x', labelsize=15 )
        ax.set_ylabel( r'$u^+_{VD}$', fontdict={'size':24} )
        ax.tick_params( axis='y', labelsize=15 )
        ax.legend( prop={'size':22,'family':'sans-serif'} ) 
        ax.set_title( r"$u^+_{VD}$ profile", size=20 )
    
    plt.savefig( figname+"_DNS" )
    plt.show()    




# ----------------------------------------------------------------------
# >>> Plot van Driest transformed u profile                      ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_u_vd :
    
    fig, ax = plt.subplots(figsize=[10,8])
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u+_vd'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
        
    ax.plot( line0.df['y+'], 
             line0.df['u+_vd'],
             line0.color, 
             label = line0.label, 
             ls    = line0.lstyle,
             linewidth = line0.width)

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
    
    figname = 'Uvd'
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
        ax.tick_params( axis='x', labelsize=15 )
        ax.set_ylabel( r'$u^+_{VD}$', fontdict={'size':24} )
        ax.tick_params( axis='y', labelsize=15 )
        ax.legend( prop={'size':22,'family':'sans-serif'} ) 
        ax.set_title( r"$u^+_{VD}$ profile", size=20 )
    
    plt.savefig( figname )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot u profile                                             ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_u :
    
    fig, ax = plt.subplots(figsize=[10,8])
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
        
    ax.plot( line0.df['y+'], 
             line0.df['u+'],
             line0.color, 
             label = line0.label, 
             ls    = line0.lstyle,
             linewidth = line0.width)

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
        ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
        ax.tick_params( axis='x', labelsize=15 )
        ax.set_ylabel( r'$u^+$', fontdict={'size':24} )
        ax.tick_params( axis='y', labelsize=15 )
        ax.legend( prop={'size':22,'family':'sans-serif'} ) 
        ax.set_title( r"$u^+$ profile", size=20 )


    plt.savefig( figname )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot Reynolds Stress Profile                              ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_RS_uu:

    fig, ax = plt.subplots( figsize=[8,8] )

    #ax.plot(dP.df['y+'],
    #        dP.df['urms+']*dP.df['urms+'],
    #        'gray',
    #        label = r'$u^\prime u^\prime$ pirozzoli',
    #        ls    = "",
    #        marker='+',
    #        \linewidth=4)

    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u`u`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
        
    ax.plot( line0.df['y+'], 
             line0.df['u`u`+'],
             line0.color, 
             label = line0.label, 
             ls    = line0.lstyle,
             linewidth = line0.width)

    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,1000] )
    ax.set_ylim( [-1,8]   )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=10,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=5,
                    width=1)
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0) )
    ax.xaxis.set_minor_locator( x_minor )

    # set spacing between major tickers.
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))

    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_uu'
    
    if pure:    
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
        ax.set_ylabel( r'$\xi u^\prime u^\prime$',
                        fontdict={'size':24} )
        ax.tick_params( axis='x', labelsize=15 )
        ax.tick_params( axis='y', labelsize=15 )
        ax.legend( prop={'size':15}, loc='upper right' ) 
        ax.set_title( "Reynolds Stress profile u`u`", size=20 )

    plt.savefig( figname )
    plt.show()


if plt_RS_vv:
    
    fig, ax = plt.subplots( figsize=[8,8] )
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['v`v`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
        
    ax.plot( line0.df['y+'], 
             line0.df['v`v`+'],
             line0.color, 
             label = line0.label, 
             ls    = line0.lstyle,
             linewidth = line0.width)
    
    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,1000] )
    ax.set_ylim( [-0.1,1.5]   )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=10,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=5,
                    width=1)
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0) )
    ax.xaxis.set_minor_locator( x_minor )
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))


    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_vv'
    
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
        ax.set_ylabel( r'$\xi v^\prime v^\prime$',
                        fontdict={'size':24} )
        ax.tick_params( axis='x', labelsize=15 )
        ax.tick_params( axis='y', labelsize=15 )
        ax.legend( prop={'size':15}, loc='upper left' ) 
        ax.set_title( "Reynolds Stress profile v`v`", size=20 )

    plt.savefig( figname )
    plt.show()
    

if plt_RS_ww:
    
    fig, ax = plt.subplots( figsize=[8,8] )
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['w`w`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
        
    ax.plot( line0.df['y+'], 
             line0.df['w`w`+'],
             line0.color, 
             label = line0.label, 
             ls    = line0.lstyle,
             linewidth = line0.width)
    
    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,1000] )
    ax.set_ylim( [-0.1,2]   )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=10,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=5,
                    width=1)
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0) )
    ax.xaxis.set_minor_locator( x_minor )
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))


    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_ww'
    
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
        ax.set_ylabel( r'$\xi w^\prime w^\prime$',
                        fontdict={'size':24} )
        ax.tick_params( axis='x', labelsize=15 )
        ax.tick_params( axis='y', labelsize=15 )
        ax.legend( prop={'size':15}, loc='upper left' ) 
        ax.set_title( "Reynolds Stress profile w`w`", size=20 )

    plt.savefig( figname )
    plt.show()


if plt_RS_uv:

    fig, ax = plt.subplots( figsize=[8,8] )
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u`v`+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
        
    ax.plot( line0.df['y+'], 
             line0.df['u`v`+'],
             line0.color, 
             label = line0.label, 
             ls    = line0.lstyle,
             linewidth = line0.width)
    
    ax.set_xscale( "symlog", linthresh=1 )


    ax.set_xlim( [1,1000] )
    ax.set_ylim( [-1,0.2] )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=10,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=5,
                    width=1)
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0) )
    ax.xaxis.set_minor_locator( x_minor )
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))


    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_uv'
    
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        figname += '_pure'
    else:
        ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
        ax.set_ylabel( r'$\xi u^\prime v^\prime$',
                        fontdict={'size':24} )
        ax.tick_params( axis='x', labelsize=15 )
        ax.tick_params( axis='y', labelsize=15 )
        ax.legend( prop={'size':15}, loc='best' ) 
        ax.set_title( "Reynolds Stress profile u`v`", size=20 )

    plt.savefig( figname )
    plt.show()
    

# ----------------------------------------------------------------------
# >>> Plot rho profile                                           ( 4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_rho :
    
    fig, ax = plt.subplots(figsize=[10,8])
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['rho'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
        
    ax.plot( line0.df['y+'], 
             line0.df['rho'],
             line0.color, 
             label = line0.label, 
             ls    = line0.lstyle,
             linewidth = line0.width)

    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
    ax.tick_params( axis='x', labelsize=15 )
    
    ax.set_ylabel( r'$rho$', fontdict={'size':24} )
    ax.tick_params( axis='y', labelsize=15 )
    
    ax.set_xlim( [1,1000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )

    ax.legend( prop={'size':20} ) 
    ax.set_title( r"$rho$ profile", size=20 )

    ax.grid()
    
    plt.savefig( "rho_shifted" )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot temperature profile                                  ( 5 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_T :
    
    fig, ax = plt.subplots(figsize=[10,8])
    


    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
    ax.tick_params( axis='x', labelsize=15 )
    
    ax.set_ylabel( r'$T/T_{\infty}$', fontdict={'size':24} )
    ax.tick_params( axis='y', labelsize=15 )
    
    ax.set_xlim( [1,1000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )

    ax.legend( prop={'size':20} ) 
#    ax.set_title( r"$u^+_{VD}$ profile", size=20 )

    ax.grid()
    
    plt.savefig( "T_shifted_new" )
    plt.show()

# ----------------------------------------------------------------------
# >>> Plot Mach_prime(fluctuating Mach) number                   ( 6 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/01/29  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Mt :
    
    fig, ax = plt.subplots(figsize=[10,8])
    


    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
    ax.tick_params( axis='x', labelsize=15 )
    
    ax.set_ylabel( r'$M_t$', fontdict={'size':24} )
    ax.tick_params( axis='y', labelsize=15 )
    
    ax.set_xlim( [1,1000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )

    ax.legend( prop={'size':20} ) 
#    ax.set_title( r"$M_t$ profile", size=20 )

    ax.grid()
    
    plt.savefig( "Mxt_profile_shifted" )
    plt.show()