# -*- coding: utf-8 -*-
'''
@File    :   plt_RSprofile.py
@Time    :   2023/02/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   comparing different RS components (for 3AF2023)
'''

import os

import numpy             as     np

import matplotlib

import matplotlib.pyplot as     plt

import matplotlib.ticker as     ticker

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

OutPath  = "/home/wencanwu/my_simulation/temp/DataPost/"

data1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/mean_result_ib_real.dat'
norm1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/statistic_average.dat'

data2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/mean_result_ib_real.dat'
norm2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/statistic_average.dat' 

data3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/mean_result_ib_real.dat'
norm3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/statistic_average.dat'

data4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/mean_result_ib_real.dat'
norm4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/statistic_average.dat'

data_piro = '/home/wencanwu/my_simulation/temp/Pirozzoli/M2_Retau_250'

data0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/x_-53.6_1storder.dat'


vd_mode  = 'bottom'


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
#d1.shift_y( 0.0 )
d1.read_norm( norm1 )
d1.get_norm()
d1.vd_transform(mode=vd_mode)
d1.get_Mt()

d2 = PlotDataframe( data2 )
d2.shift_y( 0.468 )
#d2.shift_y( 0.0 )
d2.read_norm( norm2 )
d2.get_norm()
d2.vd_transform(mode=vd_mode)
d2.get_Mt()

d3 = PlotDataframe( data3 )
d3.shift_y( 0.416 )
#d3.shift_y( 0.0 )
d3.read_norm( norm3 )
d3.get_norm()
d3.vd_transform(mode=vd_mode)
d3.get_Mt()

d4 = PlotDataframe( data4 )
d4.shift_y( 0.312 )
#d4.shift_y( 0.0 )
d4.read_norm( norm4 )
d4.get_norm()
d4.vd_transform(mode=vd_mode)
d4.get_Mt()

d0 = PlotDataframe( data0 )
d0.rho_w = d0.df['<rho>'][0]
d0.vd_transform()
d0.T_inf = 160.15
d0.df['T_nd'] = np.array(d0.df['<T>']) / d0.T_inf
d0.get_Mt()

dP = PlotDataframe( data_piro )

os.chdir( OutPath )

# ----------------------------------------------------------------------
# >>> u'u'                                                         (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/19  - created
#
# Desc
#
# ----------------------------------------------------------------------

fig, ax = plt.subplots( figsize=[8,8] )

#ax.plot(dP.df['y+'],
#        dP.df['urms+']*dP.df['urms+'],
#        'gray',
#        label = r'$u^\prime u^\prime$ pirozzoli',
#        ls    = "",
#        marker='+',
#        \linewidth=4)

ax.plot(d0.df['y_plus'],
        d0.df['<u`u`>+'],
        'gray',
        label = r'$u^\prime u^\prime$ smooth',
        ls    = "--",
        linewidth=4)
    
ax.plot(d1.df['y_s_plus'],
        d1.df['<u`u`>+'],
        'green',
        label  = r'$u^\prime u^\prime \ D/\delta_0=2.0$',
        ls     = ":",
        linewidth=4)

ax.plot(d2.df['y_s_plus'],
        d2.df['<u`u`>+'],
        'blue',
        label  = r'$u^\prime u^\prime \ D/\delta_0=1.0$',
        ls     = "-.",
        linewidth=4)

ax.plot(d3.df['y_s_plus'], 
        d3.df['<u`u`>+'],
        'black',
        label = r'$u^\prime u^\prime \ D/\delta_0=0.5$',
        ls    = (0, (3, 1, 1, 1, 1, 1)),
        linewidth=4)

ax.plot(d4.df['y_s_plus'], 
        d4.df['<u`u`>+'],
        'red',
        label = r'$u^\prime u^\prime \ D/\delta_0=0.25$',
        ls    = (0, (10, 3)),
        linewidth=4)

ax.set_xscale( "symlog", linthresh=1 )

#ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
#ax.set_ylabel( r'$\xi u^\prime u^\prime$',
#                fontdict={'size':24} )

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

#ax.tick_params( axis='x', labelsize=15 )
#ax.tick_params( axis='y', labelsize=15 )
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

#ax.legend( prop={'size':15}, loc='upper right' ) 
#ax.set_title( "Reynolds Stress profile uu", size=20 )

ax.grid(visible=True, which='both',axis='both',color='gray',
        linestyle='--',linewidth=0.2)

# Adjust the spacing around the plot to remove the white margin
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

os.chdir( OutPath )
plt.savefig( "RS_uu_pure" )
plt.show()

# ----------------------------------------------------------------------
# >>> v'v'                                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/19  - created
#
# Desc
#
# ----------------------------------------------------------------------

fig, ax = plt.subplots( figsize=[8,8] )
 
ax.plot(d0.df['y_plus'],
        d0.df['<v`v`>+'],
        'gray',
        label = r'$v^\prime v^\prime$ smooth',
        ls    = "--",
        linewidth=4)
    
ax.plot(d1.df['y_s_plus'],
        d1.df['<v`v`>+'],
        'green',
        label  = r'$v^\prime v^\prime \ D/\delta_0=2.0$',
        ls     = ":",
        linewidth=4)

ax.plot(d2.df['y_s_plus'],
        d2.df['<v`v`>+'],
        'blue',
        label  = r'$v^\prime v^\prime \ D/\delta_0=1.0$',
        ls     = "-.",
        linewidth=4)

ax.plot(d3.df['y_s_plus'], 
        d3.df['<v`v`>+'],
        'black',
        label = r'$v^\prime v^\prime \ D/\delta_0=0.5$',
        ls    = (0, (3, 1, 1, 1, 1, 1)),
        linewidth=4)

ax.plot(d4.df['y_s_plus'], 
        d4.df['<v`v`>+'],
        'red',
        label = r'$v^\prime v^\prime \ D/\delta_0=0.25$',
        ls    = (0, (10, 3)),
        linewidth=4)

ax.set_xscale( "symlog", linthresh=1 )

#ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
#ax.set_ylabel( r'$\xi v^\prime v^\prime$',
#                fontdict={'size':24} )

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

#ax.tick_params( axis='x', labelsize=15 )
#ax.tick_params( axis='y', labelsize=15 )
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

#ax.legend( prop={'size':15}, loc='upper left' ) 
#ax.set_title( "Reynolds Stress profile vv", size=20 )

ax.grid(visible=True, which='both',axis='both',color='gray',
        linestyle='--',linewidth=0.2)

# Adjust the spacing around the plot to remove the white margin
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

os.chdir( OutPath )
plt.savefig( "RS_vv_pure" )
plt.show()

# ----------------------------------------------------------------------
# >>> w'w'                                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/19  - created
#
# Desc
#
# ----------------------------------------------------------------------

fig, ax = plt.subplots( figsize=[8,8] )
 
ax.plot(d0.df['y_plus'],
        d0.df['<w`w`>+'],
        'gray',
        label = r'$w^\prime w^\prime$ smooth',
        ls    = "--",
        linewidth=4)
    
ax.plot(d1.df['y_s_plus'],
        d1.df['<w`w`>+'],
        'green',
        label  = r'$w^\prime w^\prime \ D/\delta_0=2.0$',
        ls     = ":",
        linewidth=4)

ax.plot(d2.df['y_s_plus'],
        d2.df['<w`w`>+'],
        'blue',
        label  = r'$w^\prime w^\prime \ D/\delta_0=1.0$',
        ls     = "-.",
        linewidth=4)

ax.plot(d3.df['y_s_plus'], 
        d3.df['<w`w`>+'],
        'black',
        label = r'$w^\prime w^\prime \ D/\delta_0=0.5$',
        ls    = (0, (3, 1, 1, 1, 1, 1)),
        linewidth=4)

ax.plot(d4.df['y_s_plus'], 
        d4.df['<w`w`>+'],
        'red',
        label = r'$w^\prime w^\prime \ D/\delta_0=0.25$',
        ls    = (0, (10, 3)),
        linewidth=4)

ax.set_xscale( "symlog", linthresh=1 )

#ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
#ax.set_ylabel( r'$\xi w^\prime w^\prime$',
#                fontdict={'size':24} )

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
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))

#ax.tick_params( axis='x', labelsize=15 )
#ax.tick_params( axis='y', labelsize=15 )
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

#ax.legend( prop={'size':15}, loc='upper left' ) 
#ax.set_title( "Reynolds Stress profile ww", size=20 )

ax.grid(visible=True, which='both',axis='both',color='gray',
        linestyle='--',linewidth=0.2)

# Adjust the spacing around the plot to remove the white margin
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

os.chdir( OutPath )
plt.savefig( "RS_ww_pure" )
plt.show()

# ----------------------------------------------------------------------
# >>> u'v'                                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/19  - created
#
# Desc
#
# ----------------------------------------------------------------------

fig, ax = plt.subplots( figsize=[8,8] )
 
ax.plot(d0.df['y_plus'],
        d0.df['<u`v`>+'],
        'gray',
        label = r'$u^\prime v^\prime$ smooth',
        ls    = "--",
        linewidth=4)
    
ax.plot(d1.df['y_s_plus'],
        d1.df['<u`v`>+'],
        'green',
        label  = r'$u^\prime v^\prime \ D/\delta_0=2.0$',
        ls     = ":",
        linewidth=4)

ax.plot(d2.df['y_s_plus'],
        d2.df['<u`v`>+'],
        'blue',
        label  = r'$u^\prime v^\prime \ D/\delta_0=1.0$',
        ls     = "-.",
        linewidth=4)

ax.plot(d3.df['y_s_plus'], 
        d3.df['<u`v`>+'],
        'black',
        label = r'$u^\prime v^\prime \ D/\delta_0=0.5$',
        ls    = (0, (3, 1, 1, 1, 1, 1)),
        linewidth=4)

ax.plot(d4.df['y_s_plus'], 
        d4.df['<u`v`>+'],
        'red',
        label = r'$u^\prime v^\prime \ D/\delta_0=0.25$',
        ls    = (0, (10, 3)),
        linewidth=4)

ax.set_xscale( "symlog", linthresh=1 )

#ax.set_xlabel( "$y_s^+$", fontdict={'size':24} )  
#ax.set_ylabel( r'$\xi u^\prime v^\prime$',
#                fontdict={'size':24} )

ax.set_xlim( [1,1000] )
ax.set_ylim( [-1,0.2]   )

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

#ax.tick_params( axis='x', labelsize=15 )
#ax.tick_params( axis='y', labelsize=15 )
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

#ax.legend( prop={'size':15}, loc='upper left' ) 
#ax.set_title( "Reynolds Stress profile uv", size=20 )

ax.grid(visible=True, which='both',axis='both',color='gray',
        linestyle='--',linewidth=0.2)

# Adjust the spacing around the plot to remove the white margin
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
os.chdir( OutPath )
plt.savefig( "RS_uv_pure" )
plt.show()