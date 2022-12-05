#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_profile.py
@Time    :   2022/10/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os

from   plt_tools         import *

import matplotlib.pyplot as     plt

import matplotlib

# ----------------------------------------------------------------------
# >>> Input Info                                               
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

OutPath  = "/home/wencanwu/my_simulation/temp/221125_lowRe/DataPost/"

data1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/mean_result_ib_spf.dat'
norm1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/statistic_average.dat'

data2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/mean_result_ib_spf.dat'
norm2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/statistic_average.dat' 

data3 = '/home/wencanwu/my_simulation/temp/220825_lowRe/mean_result_ib_spf.dat'
norm3 = '/home/wencanwu/my_simulation/temp/220825_lowRe/statistic_average.dat'

data4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/mean_result_ib_spf.dat'
norm4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/statistic_average.dat'

data5 = '/home/wencanwu/my_simulation/temp/221125_lowRe/mean_result_ib_spf.dat'
norm5 = '/home/wencanwu/my_simulation/temp/221125_lowRe/statistic_average.dat'

data0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/x_-68.0625.dat'

plt_u   = True
plt_RS  = True
plt_T   = False
Compare = True


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

d1 = PlotDataframe( data1 )
d1.shift_y( 0.494 )
d1.read_norm( norm1 )
d1.get_norm()

d2 = PlotDataframe( data2 )
d2.shift_y( 0.468 )
d2.read_norm( norm2 )
d2.get_norm()

d3 = PlotDataframe( data3 )
d3.shift_y( 0.416 )
d3.read_norm( norm3 )
d3.get_norm()

d4 = PlotDataframe( data4 )
d4.shift_y( 0.312 )
d4.read_norm( norm4 )
d4.get_norm()

d5 = PlotDataframe( data5 )
d5.shift_y( 0.26 )
d5.read_norm( norm5 )
d5.get_norm()

d0 = PlotDataframe( data0 )

os.chdir( OutPath )

#d0.df.to_excel("smooth.xlsx")
#d1.df.to_excel("1014.xlsx")
#d2.df.to_excel("0926.xlsx")
#d3.df.to_excel("0825.xlsx")
#d4.df.to_excel("0927.xlsx")

# ----------------------------------------------------------------------
# >>> Plot u profile                                             ( 1 )
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
    
    if Compare :
        
        ax.plot( d0.df['y_plus'], 
                 d0.df['u_plus'],
                 'gray', 
                 label = r'$smooth$', 
                 ls    = '--')

    ax.plot( d1.df['y_s_plus'], 
             d1.df['u_plus'],
             'green',   
             label = r'$D=2.0\delta_0$', 
             ls    = '-')
    
    ax.plot( d2.df['y_s_plus'], 
             d2.df['u_plus'],
             'blue',  
             label = r'$D=1.0\delta_0$', 
             ls    = '-')
    
    ax.plot( d3.df['y_s_plus'], 
             d3.df['u_plus'],
             'black', 
             label = r'$D=0.5\delta_0$', 
             ls    = '-')
    
    ax.plot( d4.df['y_s_plus'], 
             d4.df['u_plus'],
             'red',
             label  = r'$D=0.25\delta_0$', 
             ls     = '-')

    ax.plot( d5.df['y_s_plus'], 
             d5.df['u_plus'],
             'purple',
             label  = r'$D=0.125\delta_0$', 
             ls     = '-')

    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )
    ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
    ax.tick_params( axis='x', labelsize=15 )
    
    ax.set_ylabel( r'$u^+$', fontdict={'size':24} )
    ax.tick_params( axis='y', labelsize=15 )
    
    ax.set_xlim( [1,1000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base=10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )

    ax.legend( prop={'size':20} ) 
    ax.set_title( r"$u^+$ profile", size=20 )

    ax.grid()
    
    plt.savefig( "u_profile_shifted" )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot Reynolds Stress Profile                              ( 2 )
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

if plt_RS :
    
    fig, ax = plt.subplots( figsize=[10,8] )
    
    if Compare :
            
        ax.plot( d0.df['y_plus'],
                 d0.df['<u`u`>+'],
                 'blue',
                 label = r'$u^\prime u^\prime$ smooth',
                 ls    = "-")
        
        ax.plot( d0.df['y_plus'],
                 d0.df['<u`v`>+'],
                 'y',
                 label = r'$u^\prime v^\prime$ smooth',
                 ls    = "-")
        
        ax.plot( d0.df['y_plus'],
                 d0.df['<v`v`>+'],
                 'green',
                 label = r'$v^\prime v^\prime$ smooth',
                 ls    = "-")
        
        ax.plot( d0.df['y_plus'],
                 d0.df['<w`w`>+'],
                 'red',
                 label = r'$w^\prime w^\prime$ smooth',
                 ls    = "-")
    
    ax.plot( d5.df['y_s_plus'],
             d5.df['<u`u`>+'],
             'blue',
             label  = r'$u^\prime u^\prime \ D=0.125\delta_0$',
             ls     = "--")
    
    ax.plot( d5.df['y_s_plus'],
             d5.df['<u`v`>+'],
             'y',
             label  = r'$u^\prime v^\prime \ D=0.125\delta_0$',
             ls     = "--")
    
    ax.plot( d5.df['y_s_plus'], 
             d5.df['<v`v`>+'],
             'green',
             label = r'$v^\prime v^\prime \ D=0.125\delta_0$',
             ls    = "--")
    
    ax.plot( d5.df['y_s_plus'], 
             d5.df['<w`w`>+'],
             'r',
             label = r'$w^\prime w^\prime \ D=0.125\delta_0$',
             ls    = "--")

    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh=1 )
    ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
    ax.tick_params( axis='x', labelsize=15 )
    
    ax.set_xlim( [1,1000] )
    
    x_minor = matplotlib.ticker.LogLocator( 
                        base = 10.0, subs = np.arange(1.0,10.0) )
    
    ax.xaxis.set_minor_locator( x_minor )
    
    ax.set_ylabel( r'$\xi u^\prime_i u^\prime_j$',
                   fontdict={'size':24} )
    
    ax.tick_params( axis='y', labelsize=15 )

    ax.legend( prop={'size':12}, loc='upper right' ) 
    ax.set_title( "Reynolds Stress profile", size=20 )

    ax.grid()
    
    os.chdir( OutPath )
    plt.savefig( "oneeights_delta_RS" )
    plt.show()
