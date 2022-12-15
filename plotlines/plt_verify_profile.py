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


# ----------------------------------------------------------------------
# >>> Input Info                                                ( 0 )
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
