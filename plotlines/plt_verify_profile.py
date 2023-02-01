#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_varify_profile.py
@Time    :   2022/12/15 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Verify Vista can obtain the same profile as Luis
'''

import os

import numpy             as     np

import matplotlib

import matplotlib.pyplot as     plt

from   plt_tools         import PlotDataframe


# chose which ones are compared?
# mine v.s. Luis       => mode1 = True
# mine v.s. Pirozzoli  => mode2 = True

mode1 = False

mode2 = True

# ----------------------------------------------------------------------
# >>> Luis' v.s. Mine                                             ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/12/15  - created
#
# Desc
#
# ----------------------------------------------------------------------

if mode1:

    outpath = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/TBL_data_low_Re'

    data0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/x_-53.6_1storder.dat'

    dataL = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/TBL_data_low_Re/uvd_SBLI_B1.dat'


    d0 = PlotDataframe( data0 )
    d0.rho_w = d0.df['<rho>'][0]
    d0.vd_transform()
    d0.T_inf = 160.15
    d0.df['T_nd'] = np.array(d0.df['<T>']) / d0.T_inf

    dL = PlotDataframe( dataL )

    os.chdir( outpath )

    fig, ax = plt.subplots( figsize=[10,8] )
        
    ax.plot( d0.df['y_plus'], 
            d0.df['u_plus_vd'],
            'gray', 
            label = r'$smooth$', 
            ls    = '--')

    ax.plot( dL.df['y_plus'], 
            dL.df['u_plus_vd'],
            'red', 
            label = r'$from Luis$', 
            ls    = '--')


    ax.minorticks_on()

    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
    ax.tick_params( axis='x', labelsize=15 )

    ax.set_ylabel( r'$u^+_{VD}$', fontdict={'size':24} )
    ax.tick_params( axis='y', labelsize=15 )

    ax.set_xlim( [1,1000] )

    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )

    ax.xaxis.set_minor_locator( x_minor )

    ax.legend( prop={'size':20} ) 
    ax.set_title( r"$u^+_{VD}$ profile", size=20 )

    ax.grid()

    plt.savefig( "u_VD_verify_1storder" )
    plt.show()
    
# ----------------------------------------------------------------------
# >>> Mine v.s. Pirozzoli                                        ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/01/30  - created
#
# Desc
#
# ----------------------------------------------------------------------

if mode2:
    
    outpath = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/Pirozzoli'
    
    data0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/x_-53.6_1storder.dat'
#    data0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/x_-80.dat'

    data1 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/x_-80.dat'
    
    data2 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/x_-115.dat'
    
    data3 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/x_-20.dat'
    
    dataP = '/home/wencanwu/my_simulation/temp/Pirozzoli/M2_Retau_250'

    dataP2 = '/home/wencanwu/my_simulation/temp/Pirozzoli/M2_Retau_450'
    
    d0 = PlotDataframe( data0 )
    d0.rho_w = d0.df['<rho>'][0]
    d0.vd_transform()
    d0.T_inf = 160.15
    d0.df['T_nd'] = np.array(d0.df['<T>']) / d0.T_inf
    
    d1 = PlotDataframe( data1 )
    d1.rho_w = d1.df['<rho>'][0]
    d1.vd_transform()
    d1.T_inf = 160.15
    d1.df['T_nd'] = np.array(d1.df['<T>']) / d1.T_inf
    
    d2 = PlotDataframe( data2 )
    d2.rho_w = d2.df['<rho>'][0]
    d2.vd_transform()
    d2.T_inf = 160.15
    d2.df['T_nd'] = np.array(d2.df['<T>']) / d2.T_inf
    
    d3 = PlotDataframe( data3 )
    d3.rho_w = d3.df['<rho>'][0]
    d3.vd_transform()
    d3.T_inf = 160.15
    d3.df['T_nd'] = np.array(d3.df['<T>']) / d3.T_inf
    
    dP = PlotDataframe( dataP )
    dP2= PlotDataframe( dataP2 )
    
    os.chdir( outpath )
    
    fig, ax = plt.subplots( figsize=[10,8] )
    
    ax.plot( d3.df['y_plus'], 
             d3.df['u_plus_vd'],
             'yellow', 
             label = r'$smooth-20$', 
             ls    = '--')
    
    ax.plot( d0.df['y_plus'], 
             d0.df['u_plus_vd'],
             'gray', 
             label = r'$smooth-53.6$', 
             ls    = '--')
    
    ax.plot( d1.df['y_plus'], 
             d1.df['u_plus_vd'],
             'blue', 
             label = r'$smooth-80$', 
             ls    = '--')

    ax.plot( d2.df['y_plus'], 
             d2.df['u_plus_vd'],
             'green', 
             label = r'$smooth-115$', 
             ls    = '--')

    ax.plot( dP.df['y+'], 
             dP.df['u_vd+'],
             'red', 
             label = r'$Pirozzoli\_250$', 
             ls    = '-')
    
    ax.plot( dP2.df['y+'], 
             dP2.df['u_vd+'],
             'pink', 
             label = r'$Pirozzoli\_450$', 
             ls    = '-')
    
    ax.minorticks_on()

    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
    ax.tick_params( axis='x', labelsize=15 )

    ax.set_ylabel( r'$u^+_{VD}$', fontdict={'size':24} )
    ax.tick_params( axis='y', labelsize=15 )

    ax.set_xlim( [1,1000] )

    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )

    ax.xaxis.set_minor_locator( x_minor )

    ax.legend( prop={'size':20} ) 
    ax.set_title( r"$u^+_{VD}$ profile", size=20 )

    ax.grid()

    plt.savefig( "u_VD_verify_compare" )
    plt.show()