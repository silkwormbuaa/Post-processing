# -*- coding: utf-8 -*-
'''
@File    :   plt_verify_Cf.py
@Time    :   2022/12/15 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Verify my code can get variable streamwise distribution
         :   correctly.
'''
import os
import sys
import matplotlib.pyplot as     plt

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


FoldPath = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/TBL_data_low_Re'

data0_p  = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/Cf_flat_new.dat"

data_L   = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/TBL_data_low_Re/pw_SBLI_B1.dat'

delta    = 7.67277
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

d0_p = PlotDataframe( data0_p )

dL   = PlotDataframe( data_L )

d0_p.shift_x( x_imp, delta )

# ----------------------------------------------------------------------
# >>> Plot Cp                                                    ( 1 )
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
    
    fig2, ax2 = plt.subplots( figsize=[10,8.5], constrained_layout=True )

    ax2.plot( d0_p.df['x_s'], 
              d0_p.df['Cp'],
              'gray', 
              label = r'$smooth$', 
              ls    = '--')

    ax2.plot( dL.df['x_s'], 
              dL.df['Cp'],
              'red', 
              label = r'$from\ Luis$', 
              ls    = '--')

    ax2.minorticks_on()

    ax2.set_xlabel("$(x-x_{imp})/\delta_0$",fontdict={'size':40})
    ax2.tick_params(axis='x',labelsize=32)
    
    ax2.set_ylabel("$<p_w>/p_{\infty}$",fontdict={'size':40})
    ax2.tick_params(axis='y',labelsize=32)
    
    ax2.set_xlim([-10.0,12.0])
    #ax.set_ylim([0.3,0.85])

    ax2.legend(prop={'size':20})
    
    ax2.grid()

    fig2.savefig("Cp_verify.png")
    plt.show()
