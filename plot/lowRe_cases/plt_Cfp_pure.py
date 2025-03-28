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
import sys
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.plot_tools  import PlotDataframe

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


FoldPath = '/home/wencanwu/my_simulation/temp/221125_lowRe/linedata'

outpath  = '/home/wencanwu/my_simulation/temp/DataPost'

data1    = "Cf_points_1014.dat"

data2    = "Cf_points_0926.dat"

data3    = "Cf_points_0825.dat"

data4    = "Cf_points_0927.dat"

data5    = "Cf_points_1125.dat"

data0_f  = "x_cf_STBLI_Wencan.dat"

data0_p  = "Cf_flat_new.dat"

delta    = 5.2
x_imp    = 50.4

plt_Cf   = True
plt_Cp   = True

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
    
    fig1, ax1 = plt.subplots( figsize=[8,8], constrained_layout=True )
    
    
    ax1.plot( d2.df['x_s'],
              d2.df['Cf']*1000,
              'blue',
              label  = r'$D/\delta_0=1.0$',
              ls     = '-.',
              linewidth=4)
    
    ax1.plot( d3.df['x_s'],
              d3.df['Cf']*1000,
              'black',
              label  = r'$D/\delta_0=0.5$',
              ls     = (0, (3, 1, 1, 1, 1, 1)),
              linewidth=4)
    
    ax1.plot( d4.df['x_s'],
              d4.df['Cf']*1000,
              'red',
              label  = r'$D/\delta_0=0,25$', 
              ls     = (0, (10, 3)),
              linewidth=4)
    
#    ax1.plot( d5.df['x_s'],
#              d5.df['Cf']*1000,
#              'purple',
#              label  = r'$D=0.125\delta_0$', 
##              marker = "^",
#              markersize = 10)

    ax1.plot( d0_f.df['x_s'], 
              d0_f.df['Cf']*1000,
              'gray', 
              label = r'$smooth$', 
              ls    = '--',
              linewidth = 4)
    
    # dotted line are too easy to be covered, so put it at a upper layer.
    ax1.plot( d1.df['x_s'],
              d1.df['Cf']*1000,
              'green',
              label  = r'$D/\delta_0=2.0$',
              ls     = ':',
              linewidth=4)
    
    ax1.plot( [-20,12],
              [0,0],
              'black',
              ls = '--',
              linewidth=1 )
    
    
    ax1.minorticks_on()
    ax1.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=15,
                    width=2)
    ax1.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1)

    ax1.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(2.0))

    ax1.set_xlabel("$(x-x_{imp})/\delta_0$",fontdict={'size':40}) 
    ax1.tick_params(axis='x',labelsize=32)
    
    ax1.set_ylabel("$C_fx10^3$",fontdict={'size':40})
    ax1.tick_params(axis='y',labelsize=32)
    
    ax1.set_xlim([-20.0,10.0])
    ax1.set_ylim([-2.0,4.5])

#    ax1.legend(prop={'size':20},loc='best')
    

# Adjust the spacing around the plot to remove the white margin
    fig1.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax1.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    os.chdir( outpath )
    fig1.savefig("Cf_pure.png")
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
    
    fig2, ax2 = plt.subplots( figsize=[8,8], constrained_layout=True )

    
    ax2.plot( d2.df['x_s'],
              d2.df['Cp'],
              'blue',
              label  = r'$D/\delta_0=1.0$',
              ls     = '-.',
              linewidth=4)
    
    ax2.plot( d3.df['x_s'],
              d3.df['Cp'],
              'black',
              label  = r'$D/\delta_0=0.5$',
              ls     = (0, (3, 1, 1, 1, 1, 1)),
              linewidth=4)
    
    ax2.plot( d4.df['x_s'],
              d4.df['Cp'],
              'red',
              label  = r'$D/\delta_0=0.25$', 
              ls     = (0, (10, 3)),
              linewidth=4)

#    ax2.plot( d5.df['x_s'],
#              d5.df['Cp'],
#              'purple',
#              label  = r'$D=0.125\delta_0$', 
#              marker = "^",
#              markersize = 10)

    ax2.plot( d0_p.df['x_s'], 
              d0_p.df['Cp'],
              'gray', 
              label = r'$smooth$', 
              ls    = '--',
              linewidth=4)

    ax2.plot( d1.df['x_s'],
              d1.df['Cp'],
              'green',
              label  = r'$D\delta_0=2.0$',
              ls     = ':',
              linewidth=4)
    
    ax2.minorticks_on()
    ax2.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=15,
                    width=2)
    ax2.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1)

    ax2.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.4))

    ax2.set_xlabel("$(x-x_{imp})/\delta_0$",fontdict={'size':40})
    ax2.tick_params(axis='x',labelsize=32)
    
    ax2.set_ylabel("$<p_w>/p_{\infty}$",fontdict={'size':40})
    ax2.tick_params(axis='y',labelsize=32)
    
    ax2.set_xlim([-20.0,10.0])
    ax2.set_ylim([0.8,2.3])

 #   ax2.legend(prop={'size':20})
    
 # Adjust the spacing around the plot to remove the white margin
    fig2.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax2.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    fig2.savefig("Cp_pure.png")
    plt.show()



