#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_scat.py
@Time    :   2022/11/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plotting scatter with solidity
'''

import os

import numpy             as     np

import matplotlib        as     mpl

import matplotlib.pyplot as     plt

import matplotlib.ticker as     ticker

from   plt_tools         import PlotDataframe

plt.rcParams['font.family'] = "serif"

plt.rcParams.update({'font.size': 20})

OutPath = "/home/wencanwu/my_simulation/temp/DataPost"

DataFile = "statistic_compare"

xvar =  'D/δ' # or 'ESy'

Pure =  True  # if pure, without legend and label

plt_DU_vd_plus = False

plt_DU   =  False

plt_Cf   =  False

plt_vbar =  False

plt_Lsep =  True 

plt_Pmax =  True 

plt_pmax =  True 

plt_Hvor =  False 

os.chdir(OutPath)

data = PlotDataframe(DataFile)


# ----------------------------------------------------------------------
# >>> Plot DU_vd_plus                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_DU_vd_plus :
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1], 
                data.df['DU_vd+'].iloc[1],
                label=r'$\mathrm{D/\delta_0=2.0}$', 
                color='blue',
                marker='s',
                s = 100)

    ax.scatter( data.df[xvar].iloc[2], 
                data.df['DU_vd+'].iloc[2],
                label=r'$\mathrm{D/\delta_0=1.0}$', 
                color='blue',
                marker='s',
                s = 100)

    ax.scatter( data.df[xvar].iloc[3], 
                data.df['DU_vd+'].iloc[3],
                label=r'$\mathrm{D/\delta_0=0.5}$', 
                color='blue',
                marker='s',
                s = 100)

    ax.scatter( data.df[xvar].iloc[4], 
                data.df['DU_vd+'].iloc[4],
                label=r'$\mathrm{D/\delta_0=0.25}$', 
                color='blue',
                marker='s',
                s = 100)
    
    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.00,2.1] )
    ax.set_ylim( [0.80,2.4] )
        
#    ax.set_ylabel( r'$\mathrm{\Delta U_{vd}^+}$', fontdict={'size':24} )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    ax.set_axisbelow(True)
        
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
#    ax.legend(loc='best')
    
    plt.savefig("DU_vd_plus_pure.png")
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot DU                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_DU :
    
    fig, ax = plt.subplots( figsize=[10,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1], 
                data.df['DU'].iloc[1],
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s',
                s = 100)

    ax.scatter( data.df[xvar].iloc[2], 
                data.df['DU'].iloc[2],
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s',
                s = 100)

    ax.scatter( data.df[xvar].iloc[3], 
                data.df['DU'].iloc[3],
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s',
                s = 100)

    ax.scatter( data.df[xvar].iloc[4], 
                data.df['DU'].iloc[4],
                label=r'$D=0.25\delta_0$', 
                color='red',
                marker='s',
                s = 100)
    
    ax.minorticks_on()
    
#    ax.set_xscale( "log" )
    ax.set_xlabel( r'$\mathrm{ES_y}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.00,1.0] )
        
    ax.set_ylabel( r'$\mathrm{\Delta U^+}$', fontdict={'size':24} )
    
    ax.grid(True, which = 'both', ls='--')
    ax.set_axisbelow(True)
    
    ax.legend(loc='best')
    
    plt.savefig("ES_DU.png")
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot Cf                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Cf :
    
    fig, ax = plt.subplots( figsize=[10,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1], 
                data.df['Cf'].iloc[1]*1000, 
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[2], 
                data.df['Cf'].iloc[2]*1000, 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[3], 
                data.df['Cf'].iloc[3]*1000, 
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[4], 
                data.df['Cf'].iloc[4]*1000, 
                label=r'$D=0.25\delta_0$', 
                color='red',
                marker='s' ,
                s = 100)    

    ax.minorticks_on()
    
#    ax.set_xscale( "log" )
    ax.set_xlabel( r'$\mathrm{ES_y}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,1.0] )
#    ax.tick_params(axis='x',labelsize=32)
        
    ax.set_ylabel( r"$\mathrm{C_f\times 10^3}$", fontdict={'size':24} )
#    ax.tick_params(axis='y',labelsize=32)
    
    ax.grid(True, which = 'both', ls='--')
    ax.set_axisbelow(True)
    
    ax.legend(loc='best')
    
    plt.savefig("ES_Cf.png")
    plt.show()
    

# ----------------------------------------------------------------------
# >>> Plot v_max/u_infin                                         (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_vbar :
    
    fig, ax = plt.subplots( figsize=[10,10], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1], 
                data.df['v_'].iloc[1]*100, 
                label=r'$D=2.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 200)

    ax.scatter( data.df[xvar].iloc[2], 
                data.df['v_'].iloc[2]*100, 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 200)

    ax.scatter( data.df[xvar].iloc[3], 
                data.df['v_'].iloc[3]*100, 
                label=r'$D=0.5\delta_0$', 
                color='blue',
                marker='s' ,
                s = 200)

    ax.scatter( data.df[xvar].iloc[4], 
                data.df['v_'].iloc[4]*100, 
                label=r'$D=0.25\delta_0$', 
                color='blue',
                marker='s' ,
                s = 200)    

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=15,
                   width=3.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=1.5)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    
#    ax.set_xscale( "log" )
    if not Pure:
        ax.set_xlabel( r'$\mathrm{%s}$'%xvar, fontdict={'size':24} )
        ax.set_ylabel( r"$\mathrm{ v_{max}/u_{\infty}\times 10^2}$", 
                      fontdict={'size':24} )
    else:
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
#        ax.tick_params(axis='x',labelsize=32)
#        ax.tick_params(axis='y',labelsize=32)
    
    ax.set_xlim( [0.0,2.5] )
    ax.set_ylim( [1.2,3.1] )
    
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    ax.set_axisbelow(True)   # the grid lines will below the plot element.
    
    
    if not Pure:
        ax.legend(loc='best')
    if Pure:
        plt.savefig("%s_v_pure.png"%xvar[0])
    else:
        plt.savefig("%s_v_.png"%xvar[0])
        
    plt.show()
    

# ----------------------------------------------------------------------
# >>> Plot Lsep                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Lsep :
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    ax.plot( [0,2.5],
             [9.628729,9.628729],
              'black',
              ls = '--',
              linewidth=1 )

    ax.scatter( data.df[xvar].iloc[1], 
                data.df['Lsep'].iloc[1], 
                label=r'$\mathrm{D/\delta_0=2.0}$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[2], 
                data.df['Lsep'].iloc[2], 
                label=r'$\mathrm{D/\delta_0=1.0}$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[3], 
                data.df['Lsep'].iloc[3], 
                label=r'$\mathrm{D/\delta_0=0.5}$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[4], 
                data.df['Lsep'].iloc[4], 
                label=r'$\mathrm{D/\delta_0=0.25}$', 
                color='blue',
                marker='s' ,
                s = 100)    

    ax.minorticks_on()

    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.00,2.5] )
    ax.set_ylim( [7,15.0] )
        
#    ax.set_ylabel( r'$\mathrm{\Delta U_{vd}^+}$', fontdict={'size':24} )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    ax.set_axisbelow(True)
        
    plt.savefig("D_Lsep.png")
    
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
#    ax.legend(loc='best')
    
    plt.savefig("D_Lsep_pure.png")
    
    plt.show()
    

# ----------------------------------------------------------------------
# >>> Plot Pmax (maximum pressure after SW)                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Pmax :
    
    fig, ax = plt.subplots( figsize=[10,8], 
                            constrained_layout=True )

    ax.plot( [0,2.5],
             [2.2821733,2.2821733],
              'black',
              ls = '--',
              linewidth=1 )
    
    ax.scatter( data.df[xvar].iloc[1], 
                data.df['Pmax'].iloc[1], 
                label=r'$D/\delta_0=2.0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[2], 
                data.df['Pmax'].iloc[2], 
                label=r'$D/\delta_0=1.0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[3], 
                data.df['Pmax'].iloc[3], 
                label=r'$D/\delta_0=0.5$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[4], 
                data.df['Pmax'].iloc[4], 
                label=r'$D/\delta_0=0.25$', 
                color='blue',
                marker='s' ,
                s = 100)    

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.050))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.00,2.5] )
    ax.set_ylim( [2.15,2.3] )
        
#    ax.set_ylabel( r'$\mathrm{\Delta U_{vd}^+}$', fontdict={'size':24} )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    ax.set_axisbelow(True)
        
    plt.savefig("D_Pmax.png")
    
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
#    ax.legend(loc='best')
    
    plt.savefig("D_Pmax_pure.png")
    plt.show()
    

# ----------------------------------------------------------------------
# >>> Plot pmax ( maximum pressure fluctuation )                 (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_pmax :
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    ax.plot( [0,2.5],
             [0.083018,0.083018],
              'black',
              ls = '--',
              linewidth=1 )


    ax.scatter( data.df[xvar].iloc[1], 
                data.df['p`max'].iloc[1], 
                label=r'$D/\delta_0=2.0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[2], 
                data.df['p`max'].iloc[2], 
                label=r'$D/\delta_0=1.0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[3], 
                data.df['p`max'].iloc[3], 
                label=r'$D/\delta_0=0.5$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[4], 
                data.df['p`max'].iloc[4], 
                label=r'$D/\delta_0=0.25$', 
                color='blue',
                marker='s' ,
                s = 100)    

    ax.set_xlabel( r'$\mathrm{ ES_y }$', fontdict={'size':24} )
    ax.set_ylabel( r"$\mathrm{ p'_{max}/p_{\infty} }$", fontdict={'size':24} )

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.010))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.00,2.5] )
    ax.set_ylim( [0.07,0.1] )
        
#    ax.set_ylabel( r'$\mathrm{\Delta U_{vd}^+}$', fontdict={'size':24} )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    ax.set_axisbelow(True)
    
    plt.savefig("D_pfluc.png")
        
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
#    ax.legend(loc='best')
    
    plt.savefig("D_pfluc_pure.png")
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot Hvor ( height or size of secondary vortices )          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Hvor :
    
    fig, ax = plt.subplots( figsize=[10,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1], 
                data.df['Hvor'].iloc[1], 
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[2], 
                data.df['Hvor'].iloc[2], 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[3], 
                data.df['Hvor'].iloc[3], 
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s' ,
                s = 100)

    ax.scatter( data.df[xvar].iloc[4], 
                data.df['Hvor'].iloc[4], 
                label=r'$D=0.25\delta_0$', 
                color='red',
                marker='s' ,
                s = 100)    

    ax.minorticks_on()
    
#    ax.set_xscale( "log" )
    ax.set_xlabel( r'$\mathrm{ ES_y }$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,1.0] )
#    ax.tick_params(axis='x',labelsize=32)
        
    ax.set_ylabel( r"$\mathrm{ H/\delta_0 }$", fontdict={'size':24} )
#    ax.tick_params(axis='y',labelsize=32)
    
    ax.grid(True, which = 'both', ls='--')
    ax.set_axisbelow(True)
    
    ax.legend(loc='best')
    
    plt.savefig("ES_Hvor.png")
    plt.show()    