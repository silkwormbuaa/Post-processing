#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_statistics.py
@Time    :   2024/04/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.plot_tools   import PlotDataframe

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size'] = 40

# =============================================================================

Datapath = "/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe_ridge_height"
outpath = "/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe_ridge_height/statistics"


DataFile = "statistic_compare"

xvar =  'H/δ' # or 'ESy'

pure =  False  # if pure, without legend and label

fmt = '.pdf'  # '.pdf'

show = False

plt_DU_vd_plus =  True      # roughness function based on vd transformed velocity
plt_Cf         =  False      # skin friction coefficient
plt_vbar       =  False      # normalized vertical velocity
plt_Lsep       =  False      # length of separation
plt_Asep       =  False      # area of separation
plt_Pmax       =  False      # maximum wall pressure
plt_pmax       =  False      # maximum pressure fluctuation
plt_Hvor       =  False      # height of vortex
plt_Hson       =  False      # height of sonic line
plt_pt         =  False      # total pressure change
plt_bubble     =  False      # bubble size
plt_bubble_dev =  False      # bubble size deviation


os.chdir(Datapath)

data = PlotDataframe(DataFile)

os.chdir( outpath )


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

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['DU_vd+'].iloc[1:],
                color='blue',
                marker='s',
                s = 200 )
    
    ax.plot( [0.0,0.30],[0.0,0.0],
             color='gray',
             linestyle='--',
             linewidth=4.0)
    
    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0, 0.25] )
    ax.set_ylim( [-0.5,3.5] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
        
    figname = "DU_vd_plus"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r'$\Delta \langle u \rangle_{vD}^+$', labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    ax.set_axisbelow(False)  # set axis order at bottom

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )
    
    if show: plt.show()


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
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['Cf'].iloc[1:]*1000, 
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)
    
    Cf_smooth = data.df['Cf'].iloc[0]*1000
    ax.plot( [0.0,0.30],[Cf_smooth,Cf_smooth],
             color='gray',
             linestyle='--',
             linewidth=4.0,
             zorder=0)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [2.8,4.2] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
        
    figname = "Cf"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r'$C_f\cdot 1000$', labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )

    if show: plt.show()    

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
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['v_'].iloc[1:]*100, 
                color='blue',
                marker='s' ,
                s = 200 )

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [0.0,3.5] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
    figname = "v_max"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r'$\langle v  \rangle _{max} /u_{\infty} \cdot 100$', 
                       labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )
    
    if show: plt.show()    

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

    Lsep_smooth = data.df['Lsep'].iloc[0]
    ax.plot( [0.0,0.30],
             [Lsep_smooth,Lsep_smooth],
             'gray',
             ls = '--',
             linewidth=4.0,
             zorder=0)

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['Lsep'].iloc[1:], 
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [8.0,18.0] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
    figname = "Lsep"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r'$L_{sep}/\delta_0$', labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )

    if show: plt.show()    

# ----------------------------------------------------------------------
# >>> Plot Separation Area                                       (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/22  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Asep :
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    Asep_smooth = data.df['Asep'].iloc[0]
    ax.plot( [0.0,0.30],
             [Asep_smooth,Asep_smooth],
             'gray',
             ls = '--',
             linewidth=4.0,
             zorder=0)

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['Asep'].iloc[1:], 
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [0.8,2.0] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
    figname = "Asep"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r'$A_{sep}/A_{sep_{smooth}}$', labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )

    if show: plt.show()


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
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    Pmax_smooth = data.df['Pmax'].iloc[0]
    ax.plot( [0,2.5],
             [Pmax_smooth,Pmax_smooth],
             'gray',
             ls = '--',
             linewidth=4.0,
             zorder=0)
    
    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['Pmax'].iloc[1:], 
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [2.05,2.3] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
        
    figname = "Pmax"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r'$\langle p \rangle_{max}/p_{\infty}$', 
                       labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )

    if show: plt.show()

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

    pfluc_max_smooth = data.df['p`max'].iloc[0]
    ax.plot( [0,2.5],
             [pfluc_max_smooth,pfluc_max_smooth],
              'gray',
              ls = '--',
              linewidth=4,
              zorder=0)

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['p`max'].iloc[1:], 
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.01))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [0.05,0.09] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
    figname = "pfluc_max"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r"$\sqrt{\langle p^{'}p^{'} \rangle_{max}}/p_{\infty}$", 
                       labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )

    if show: plt.show()

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
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    # use local boundary layer thickness as reference
    
    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['Hvor'].iloc[1:]*5.2/6.02925,  
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [0.15,0.35] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
    figname = "Hvor"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r"$H_{vortex}/\delta_{99}$", 
                       labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )

    if show: plt.show()    

# ----------------------------------------------------------------------
# >>> plot Hson (height of sonic line)                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/24  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Hson :
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['H_sonic'].iloc[1:] *5.2/6.02925, 
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)

    Hson_smooth = data.df['H_sonic'].iloc[0] *5.2/6.02925
    ax.plot( [0,2.5],
             [Hson_smooth,Hson_smooth],
              'gray',
              ls = '--',
              linewidth=4,
              zorder=0)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [0.05,0.20] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
    figname = "Hson"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r"$H_{sonic}/\delta_{99}$", 
                       labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )

    if show: plt.show()
    

# ----------------------------------------------------------------------
# >>> Plot total pressure recovery                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/27  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_pt:

    fig, ax = plt.subplots( figsize=[8.8,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['pt_recovery'].iloc[1:], 
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)
    
    ax.plot( [0.0,3.0],[0.9391,0.9391],
             color='gray',
             linestyle='--',
             linewidth=4.0)
    
    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.005))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [0.935,0.955] )
    
#    ax.grid(visible=True, which='major',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
    
    figname = "pt_recovery"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r"$p_t/p_{t_0}$", 
                       labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+fmt )

    if show: plt.show()
    

# ----------------------------------------------------------------------
# >>> bubble size                                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/05/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_bubble:
    
    fig,ax = plt.subplots( figsize=[8.8,8], 
                           constrained_layout=True )
    
    bubble_smooth = data.df['bubble_stat'].iloc[0]
    
    ax.scatter( data.df[xvar].iloc[1:],
                data.df['bubble_stat'].iloc[1:]/bubble_smooth,
                color='blue',
                marker='s',
                s=200,
                zorder=10)
    
    ax.plot( [0.0,3.0],[1.0,1.0],
             color='gray',
             linestyle='--',
             linewidth=4.0)
    
    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [0.0,5.0] )
    
    figname = 'bubble_stat'
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r"$\overline{V}/\overline{V}_{smooth}$", 
                       labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )
        
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(3)
        
        plt.savefig( figname+fmt )
        
        if show: plt.show()
        
        plt.close()
        

# ----------------------------------------------------------------------
# >>> bubble size deviation                                     (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/05/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_bubble_dev:
    
    fig,ax = plt.subplots( figsize=[8.8,8], 
                           constrained_layout=True )
    
    bubble_smooth = data.df['bubble_size'].iloc[0]
    
    ax.scatter( data.df[xvar].iloc[1:],
                data.df['bubble_dev'].iloc[1:]/data.df['bubble_size'].iloc[1:],
                color='blue',
                marker='s',
                s=200,
                zorder=10)
    
    ax.plot( [0.0,3.0],[data.df['bubble_dev'].iloc[0]/bubble_smooth,data.df['bubble_dev'].iloc[0]/bubble_smooth],
             color='gray',
             linestyle='--',
             linewidth=4.0)
    
    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=3.0)
    
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=2.0)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
    
    ax.set_xlim( [0.0,0.25] )
    ax.set_ylim( [0.04,0.12] )
    
    figname = 'bubble_dev_norm'
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    else:
        ax.set_xlabel( r"$H/\delta_0$", labelpad=0 )
        ax.set_ylabel( r"$\sigma_{V/\overline{V}}$", 
                       labelpad=0 )
        ax.tick_params( axis='x', pad=15 )
        ax.tick_params( axis='y', pad=10 )
        
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(3)
        
        plt.savefig( figname+fmt )
        
        if show: plt.show()
        
        plt.close()