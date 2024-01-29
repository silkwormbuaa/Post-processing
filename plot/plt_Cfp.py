# -*- coding: utf-8 -*-
'''
@File    :   plt_C.py
@Time    :   2022/11/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Script for plotting Cf and Cp
'''


import os

import pickle

import numpy             as     np

import matplotlib.pyplot as     plt

from   plt_tools         import PlotDataframe


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


FoldPath = '/media/wencanwu/Seagate Expansion Drive/temp/221221/linedata'

outpath  = '/media/wencanwu/Seagate Expansion Drive/temp/221221/linedata'


data1    = "Cf_points_1014.dat"

data2    = "Cf_points_0926.dat"

data3    = "Cf_points_0825.dat"

data4    = "Cf_points_0927.dat"

data5    = "Cf_points_1221.dat"

data0_f  = "x_cf_STBLI_Wencan.dat"

data0_p  = "Cf_flat_new.dat"


delta    = 5.2
x_imp    = 50.4

plt_Cf   = True
plt_Cp   = True

df1file = '/media/wencanwu/Seagate Expansion Drive/temp/221014/results/streamwise_vars.pkl'
df2file = '/media/wencanwu/Seagate Expansion Drive/temp/220926/results/streamwise_vars.pkl'
df3file = '/media/wencanwu/Seagate Expansion Drive/temp/220825/results/streamwise_vars.pkl'
df4file = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/streamwise_vars.pkl'
df5file = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results/streamwise_vars.pkl'
with open(df1file,'rb') as f:  df1 = pickle.load( f )
with open(df2file,'rb') as f:  df2 = pickle.load( f )
with open(df3file,'rb') as f:  df3 = pickle.load( f )
with open(df4file,'rb') as f:  df4 = pickle.load( f )
with open(df5file,'rb') as f:  df5 = pickle.load( f )

# ----------------------------------------------------------------------
# >>> Initialize data                                            ( 0 )
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

os.chdir( FoldPath )

d1   = PlotDataframe( data1 )
d2   = PlotDataframe( data2 )
d3   = PlotDataframe( data3 )
d4   = PlotDataframe( data4 )
d5   = PlotDataframe( data5 )
d0_f = PlotDataframe( data0_f )
d0_p = PlotDataframe( data0_p )

d1.shift_x( x_imp, delta )
d2.shift_x( x_imp, delta )
d3.shift_x( x_imp, delta )
d4.shift_x( x_imp, delta )
d5.shift_x( x_imp, delta )
d0_f.shift_x( x_imp, delta )
d0_p.shift_x( x_imp, delta )

# ----------------------------------------------------------------------
# >>> Plot Cf                                                    ( 1 )
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

if plt_Cf is True:
    
    fig1, ax1 = plt.subplots( figsize=[10,8.5], constrained_layout=True )
    
    ax1.plot( d5.df['x_s'],
              d5.df['Cf']*1000,
              'purple',
              label  = r'$D/\delta_0=0.125$', 
#              marker = "^",
              markersize = 10,
              linewidth=4)
    
    ax1.plot( d4.df['x_s'],
              d4.df['Cf']*1000,
              'red',
              label  = r'$D/\delta_0=0.25$', 
              ls     = (0, (7, 3)),
              linewidth=4)
    
    ax1.plot( d3.df['x_s'],
              d3.df['Cf']*1000,
              'black',
              label  = r'$D/\delta_0=0.5$',
              ls     = (0, (3, 1, 1, 1, 1, 1)),
              linewidth=4)
    
    ax1.plot( d2.df['x_s'],
              d2.df['Cf']*1000,
              'blue',
              label  = r'$D/\delta_0=1.0$',
              ls     = '-.',
              linewidth=4)

    ax1.plot( d1.df['x_s'],
              d1.df['Cf']*1000,
              'green',
              label  = r'$D/\delta_0=2.0$',
              ls     = ':',
              linewidth=4)
    
    ax1.plot( d0_f.df['x_s'], 
              d0_f.df['Cf']*1000,
              'gray', 
              label = r'$smooth$', 
              ls    = '--',
              linewidth = 4)
    
   
    ax1.plot( [-20,12],
              [0,0],
              'black',
              ls = '--',
              linewidth=1 )

    ax1.plot( df1['x'], df1['Cf'], 'green',ls='--')
    ax1.plot( df2['x'], df2['Cf'], 'blue', ls='--')
    ax1.plot( df3['x'], df3['Cf'], 'black',ls='--')
    ax1.plot( df4['x'], df4['Cf'], 'red',ls='--')    
    ax1.plot( df5['x'], df5['Cf'], 'purple',ls='--')
    
    ax1.minorticks_on()

    ax1.set_xlabel("$(x-x_{imp})/\delta_0$",fontdict={'size':40}) 
    ax1.tick_params(axis='x',labelsize=32)
    
    ax1.set_ylabel("$C_fx10^3$",fontdict={'size':40})
    ax1.tick_params(axis='y',labelsize=32)
    
    ax1.set_xlim([-20.0,10.0])
    ax1.set_ylim([-2.0,4.5])

    ax1.legend(prop={'size':20},loc='best')
    
    ax1.grid()

    fig1.savefig("Cf_compare.png")
    plt.show()

    
# ----------------------------------------------------------------------
# >>> Plot Cp                                                     ( 2 )
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

if plt_Cp is True:
    
    fig2, ax2 = plt.subplots( figsize=[12,9], constrained_layout=True )

    ax2.plot( d5.df['x_s'],
              d5.df['Cp'],
              'purple',
              label  = r'$D/\delta_0=0.125$', 
#              marker = "^",
#              markersize = 10,
              linewidth=4)
    
#    ax2.plot( [-10.6786859,-10.6786859],
#              [0.8,2.4],
#              'purple',
#              linewidth=1)
    
    ax2.plot( d4.df['x_s'],
              d4.df['Cp'],
              'red',
              label  = r'$D/\delta_0=0.25$', 
              ls     = (0, (7, 3)),
              linewidth=4)

#    ax2.plot( [-10.6574429,-10.6574429],
#              [0.8,2.4],
#              'red',
#              linewidth=1)
    
    ax2.plot( d3.df['x_s'],
              d3.df['Cp'],
              'black',
              label  = r'$D/\delta_0=0.5$',
              ls     = (0, (3, 1, 1, 1, 1, 1)),
              linewidth=4)

#    ax2.plot( [-9.179693795,-9.179693795],
#              [0.8,2.4],
#              'black',
#              linewidth=1)

    ax2.plot( d2.df['x_s'],
              d2.df['Cp'],
              'blue',
              label  = r'$D/\delta_0=1.0$',
              ls     = '-.',
              linewidth=4)

#    ax2.plot( [-8.316817364,-8.316817364],
#              [0.8,2.4],
#              'blue',
#              linewidth=1)
    
    ax2.plot( d1.df['x_s'],
              d1.df['Cp'],
              'green',
              label  = r'$D/\delta_0=2.0$',
              ls     = ':',
              linewidth=4)

#    ax2.plot( [-8.405281852,-8.405281852],
#              [0.8,2.4],
#              'green',
#              linewidth=1)

    ax2.plot( d0_p.df['x_s'], 
              d0_p.df['Cp'],
              'gray', 
              label = r'$smooth$', 
              ls    = '--',
              linewidth=4)
    
    ax2.plot( df4['x'], df4['Cp'], 'yellow', ls ='--')

#    ax2.plot( [-8.56077913,-8.56077913],
#              [0.8,2.4],
#              'gray',
#              linewidth=1)

    ax2.plot( df1['x'], df1['Cp'], 'green',ls='--')
    ax2.plot( df2['x'], df2['Cp'], 'blue', ls='--')
    ax2.plot( df3['x'], df3['Cp'], 'black',ls='--')
    ax2.plot( df4['x'], df4['Cp'], 'red',ls='--')    
    ax2.plot( df5['x'], df5['Cp'], 'purple',ls='--')

    ax2.minorticks_on()

    ax2.set_xlabel("$(x-x_{imp})/\delta_0$",fontdict={'size':40})
    ax2.tick_params(axis='x',labelsize=32)
    
    ax2.set_ylabel("$<p_w>/p_{\infty}$",fontdict={'size':40})
    ax2.tick_params(axis='y',labelsize=32)
    
    ax2.set_xlim([-20.0,10.0])
    ax2.set_ylim([0.8,2.4])

    ax2.legend(prop={'size':20})
    
    ax2.grid()

    fig2.savefig("Cp_compare.png")
    plt.show()



