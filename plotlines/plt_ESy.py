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

from   plt_tools         import *

plt.rcParams['font.family'] = "serif"

plt.rcParams.update({'font.size': 20})

OutPath = "/home/wencanwu/my_simulation/temp/221014_lowRe/DataPost"

DataFile = "statistic_compare"

plt_DU   =  True

plt_Cf   =  True

plt_vbar =  True

plt_Lsep =  True 

plt_Pmax =  True 

plt_pmax =  True 

plt_Hvor =  True 

os.chdir(OutPath)

data = PlotDataframe(DataFile)


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

    ax.scatter( data.df['ESy'].iloc[1], 
                data.df['DU'].iloc[1],
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s',
                s = 100)

    ax.scatter( data.df['ESy'].iloc[2], 
                data.df['DU'].iloc[2],
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s',
                s = 100)

    ax.scatter( data.df['ESy'].iloc[3], 
                data.df['DU'].iloc[3],
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s',
                s = 100)

    ax.scatter( data.df['ESy'].iloc[4], 
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

    ax.scatter( data.df['ESy'].iloc[1], 
                data.df['Cf'].iloc[1]*1000, 
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[2], 
                data.df['Cf'].iloc[2]*1000, 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[3], 
                data.df['Cf'].iloc[3]*1000, 
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[4], 
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
    
    fig, ax = plt.subplots( figsize=[10,8], 
                            constrained_layout=True )

    ax.scatter( data.df['ESy'].iloc[1], 
                data.df['v_'].iloc[1]*100, 
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[2], 
                data.df['v_'].iloc[2]*100, 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[3], 
                data.df['v_'].iloc[3]*100, 
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[4], 
                data.df['v_'].iloc[4]*100, 
                label=r'$D=0.25\delta_0$', 
                color='red',
                marker='s' ,
                s = 100)    

    ax.minorticks_on()
    
#    ax.set_xscale( "log" )
    ax.set_xlabel( r'$\mathrm{ES_y}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,1.0] )
#    ax.tick_params(axis='x',labelsize=32)
        
    ax.set_ylabel( r"$\mathrm{ v_{max}/u_{\infty}\times 10^2}$", fontdict={'size':24} )
#    ax.tick_params(axis='y',labelsize=32)
    
    ax.grid(True, which = 'both', ls='--')
    ax.set_axisbelow(True)
    
    ax.legend(loc='best')
    
    plt.savefig("ES_v_.png")
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
    
    fig, ax = plt.subplots( figsize=[10,8], 
                            constrained_layout=True )

    ax.scatter( data.df['ESy'].iloc[1], 
                data.df['Lsep'].iloc[1], 
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[2], 
                data.df['Lsep'].iloc[2], 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[3], 
                data.df['Lsep'].iloc[3], 
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[4], 
                data.df['Lsep'].iloc[4], 
                label=r'$D=0.25\delta_0$', 
                color='red',
                marker='s' ,
                s = 100)    

    ax.minorticks_on()
    
#    ax.set_xscale( "log" )
    ax.set_xlabel( r'$\mathrm{ ES_y }$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,1.0] )
#    ax.tick_params(axis='x',labelsize=32)
        
    ax.set_ylabel( r"$\mathrm{ L_{sep}/\delta_0 }$", fontdict={'size':24} )
#    ax.tick_params(axis='y',labelsize=32)
    
    ax.grid(True, which = 'both', ls='--')
    ax.set_axisbelow(True)
    
    ax.legend(loc='best')
    
    plt.savefig("ES_Lsep.png")
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

    ax.scatter( data.df['ESy'].iloc[1], 
                data.df['Pmax'].iloc[1], 
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[2], 
                data.df['Pmax'].iloc[2], 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[3], 
                data.df['Pmax'].iloc[3], 
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[4], 
                data.df['Pmax'].iloc[4], 
                label=r'$D=0.25\delta_0$', 
                color='red',
                marker='s' ,
                s = 100)    

    ax.minorticks_on()
    
#    ax.set_xscale( "log" )
    ax.set_xlabel( r'$\mathrm{ ES_y }$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,1.0] )
#    ax.tick_params(axis='x',labelsize=32)
        
    ax.set_ylabel( r"$\mathrm{ p_{max}/p_{\infty} }$", fontdict={'size':24} )
#    ax.tick_params(axis='y',labelsize=32)
    
    ax.grid(True, which = 'both', ls='--')
    ax.set_axisbelow(True)
    
    ax.legend(loc='best')
    
    plt.savefig("ES_Pmax.png")
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
    
    fig, ax = plt.subplots( figsize=[10,8], 
                            constrained_layout=True )

    ax.scatter( data.df['ESy'].iloc[1], 
                data.df['p`max'].iloc[1], 
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[2], 
                data.df['p`max'].iloc[2], 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[3], 
                data.df['p`max'].iloc[3], 
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[4], 
                data.df['p`max'].iloc[4], 
                label=r'$D=0.25\delta_0$', 
                color='red',
                marker='s' ,
                s = 100)    

    ax.minorticks_on()
    
#    ax.set_xscale( "log" )
    ax.set_xlabel( r'$\mathrm{ ES_y }$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,1.0] )
#    ax.tick_params(axis='x',labelsize=32)
        
    ax.set_ylabel( r"$\mathrm{ p'_{max}/p_{\infty} }$", fontdict={'size':24} )
#    ax.tick_params(axis='y',labelsize=32)
    
    ax.grid(True, which = 'both', ls='--')
    ax.set_axisbelow(True)
    
    ax.legend(loc='best')
    
    plt.savefig("ES_pmax.png")
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

    ax.scatter( data.df['ESy'].iloc[1], 
                data.df['Hvor'].iloc[1], 
                label=r'$D=2.0\delta_0$', 
                color='green',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[2], 
                data.df['Hvor'].iloc[2], 
                label=r'$D=1.0\delta_0$', 
                color='blue',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[3], 
                data.df['Hvor'].iloc[3], 
                label=r'$D=0.5\delta_0$', 
                color='black',
                marker='s' ,
                s = 100)

    ax.scatter( data.df['ESy'].iloc[4], 
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