#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_uprofile.py
@Time    :   2023/09/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

import numpy             as     np

import matplotlib

import matplotlib.pyplot as     plt

from   plt_tools         import PlotDataframe

from   matplotlib.font_manager import FontProperties

data1 = '/home/wencanwu/my_simulation/temp/220927_lowRe/mean_result_ib_real.dat'

data2 = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/output.txt'

norm = '/home/wencanwu/my_simulation/temp/220927_lowRe/statistic_average.dat'

vd_mode = 'bottom'

d1 = PlotDataframe( data1 )
d1.shift_y( 0.312 )
d1.read_norm( norm )
d1.get_norm()
d1.vd_transform(mode=vd_mode)
d1.get_Mt()


d2 = PlotDataframe( data2 )
d2.shift_y( 0.312 )
d2.read_norm( norm )
d2.get_norm()
d2.vd_transform(mode=vd_mode)
d2.get_Mt()

fig, ax = plt.subplots(figsize=[10,8])



ax.plot( d1.df['y_s_plus'], 
            d1.df['u_plus'],
            'green',   
            label = r'$old D=0.25\delta_0$', 
            ls    = '-')

ax.plot( d2.df['y_s_plus'], 
            d2.df['u_plus'],
            'blue',  
            label = r'$new D=0.25\delta_0$', 
            ls    = '-')




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

#plt.savefig( "u_profile_shifted" )
plt.show()

