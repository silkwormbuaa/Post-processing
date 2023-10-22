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

import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker

from   plt_tools          import PlotDataframe

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{bm}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams.update({'font.size': 20})


# =============================================================================

Datapath = "/home/wencanwu/my_simulation/temp/DataPost/"

DataFile = "statistic_compare"

xvar =  'D/δ' # or 'ESy'

pure =  True  # if pure, without legend and label

plt_DU_vd_plus = True  # roughness function based on vd transformed velocity

plt_Cf   =  True     # skin friction coefficient

plt_vbar =  True     # normalized vertical velocity

plt_Lsep =  True     # length of separation

plt_Asep =  True      # area of separation

plt_Pmax =  True     # maximum wall pressure

plt_pmax =  True     # maximum pressure fluctuation

plt_Hvor =  True     # height of vortex

os.chdir(Datapath)

data = PlotDataframe(DataFile)

os.chdir('./statistics')


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
    
    ax.plot( [0.1,3.0],[0.0,0.0],
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
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,2.1] )
    ax.set_ylim( [-0.5,2.5] )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
        
    figname = "DU_vd_plus"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$D/\delta_0$", fontsize=30, labelpad=0 )
        ax.set_ylabel( r'$\Delta U_{vd}^+$', fontsize=30, labelpad=0 )
        ax.tick_params( axis='x', labelsize=30, pad=15 )
        ax.tick_params( axis='y', labelsize=30, pad=10 )

    ax.set_axisbelow(False)  # set axis order at bottom

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+'.pdf' )
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
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['Cf'].iloc[1:]*1000, 
                color='blue',
                marker='s' ,
                s = 200,
                zorder=10)
    
    Cf_smooth = data.df['Cf'].iloc[0]*1000
    ax.plot( [0.0,3.0],[Cf_smooth,Cf_smooth],
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
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,2.1] )
    ax.set_ylim( [2.7,3.5] )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
        
    figname = "Cf"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$D/\delta_0$", fontsize=30, labelpad=0 )
        ax.set_ylabel( r'$C_f\cdot 1000$', fontsize=30, labelpad=0 )
        ax.tick_params( axis='x', labelsize=30, pad=15 )
        ax.tick_params( axis='y', labelsize=30, pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+'.pdf' )
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
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,2.1] )
    ax.set_ylim( [0.0,3.5] )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    figname = "v_max"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$D/\delta_0$", fontsize=30, labelpad=0 )
        ax.set_ylabel( r'$\langle v_{max} \rangle /u_{ref} \cdot 100$', 
                       fontsize=30, labelpad=0 )
        ax.tick_params( axis='x', labelsize=30, pad=15 )
        ax.tick_params( axis='y', labelsize=30, pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+'.pdf' )
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
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    Lsep_smooth = data.df['Lsep'].iloc[0]
    ax.plot( [0.0,2.5],
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
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,2.1] )
    ax.set_ylim( [8.0,14.0] )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    figname = "Lsep"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$D/\delta_0$", fontsize=30, labelpad=0 )
        ax.set_ylabel( r'$L_{sep}/\delta_0$', fontsize=30, labelpad=0 )
        ax.tick_params( axis='x', labelsize=30, pad=15 )
        ax.tick_params( axis='y', labelsize=30, pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+'.pdf' )
    plt.show()
    

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
    ax.plot( [0.0,2.5],
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
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,2.1] )
    ax.set_ylim( [0.8,1.6] )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    figname = "Asep"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$D/\delta_0$", fontsize=30, labelpad=0 )
        ax.set_ylabel( r'$A_{sep}/A_{sep_{smooth}}$', fontsize=30, labelpad=0 )
        ax.tick_params( axis='x', labelsize=30, pad=15 )
        ax.tick_params( axis='y', labelsize=30, pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+'.pdf' )
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
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,2.1] )
    ax.set_ylim( [2.15,2.3] )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
        
    figname = "Pmax"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$D/\delta_0$", fontsize=30, labelpad=0 )
        ax.set_ylabel( r'$\langle p \rangle_{max}/p_{ref}$', 
                       fontsize=30, labelpad=0 )
        ax.tick_params( axis='x', labelsize=30, pad=15 )
        ax.tick_params( axis='y', labelsize=30, pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+'.pdf' )
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
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.01))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,2.1] )
    ax.set_ylim( [0.07,0.10] )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    figname = "pfluc_max"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$D/\delta_0$", fontsize=30, labelpad=0 )
        ax.set_ylabel( r"$\sqrt{\langle p^{'}p^{'} \rangle_{max}}/p_{ref}$", 
                       fontsize=30, labelpad=0 )
        ax.tick_params( axis='x', labelsize=30, pad=15 )
        ax.tick_params( axis='y', labelsize=30, pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+'.pdf' )
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
    
    fig, ax = plt.subplots( figsize=[8,8], 
                            constrained_layout=True )

    ax.scatter( data.df[xvar].iloc[1:], 
                data.df['Hvor'].iloc[1:], 
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
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    
#    ax.set_xscale( "log" )
#    ax.set_xlabel( r'$\mathrm{D/delta_0}$', fontdict={'size':24} )
    
    ax.set_xlim( [0.0,2.1] )
    ax.set_ylim( [0.0,0.8] )
    
    ax.grid(visible=True, which='major',axis='both',color='gray',
            linestyle='--',linewidth=0.5)
    
    figname = "Hvor"
    
    if pure:
        figname = figname + '_pure'
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
    else:
        ax.set_xlabel( r"$D/\delta_0$", fontsize=30, labelpad=0 )
        ax.set_ylabel( r"$H_{vortice}/\delta_{0}$", 
                       fontsize=30, labelpad=0 )
        ax.tick_params( axis='x', labelsize=30, pad=15 )
        ax.tick_params( axis='y', labelsize=30, pad=10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname+'.pdf' )
    plt.show()