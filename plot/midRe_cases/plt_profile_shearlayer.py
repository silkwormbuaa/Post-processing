#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_profile_shearlayer.py
@Time    :   2025/05/15 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plot the profiles in outer scale, mainly to compare the difference in shear layer.
'''

import os
import sys
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line         import ProfileData
from   vista.directories  import create_folder

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 40

# =============================================================================

plt_u              = False
plt_RS_uu          = False
plt_RS_vv          = False
plt_RS_ww          = False
plt_RS_uv          = False
plt_T              = False
plt_Tt             = False
plt_pt             = False
plt_tke            = True
plt_p_fluc         = False

location = 'shearlayer'   # 'upstream' or 'incip' 

fmt = '.png'

# =============================================================================

outPath  = f'/home/wencan/temp/DataPost/midRe/profile_{location}/'

data0 = f'/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/profile_{location}'
data1 = f'/media/wencan/Expansion/temp/smooth_mid/postprocess/statistics/profile_{location}'
data2 = f'/media/wencan/Expansion/temp/220927/postprocess/statistics/profile_{location}'
data3 = f'/media/wencan/Expansion/temp/241030/postprocess/statistics/profile_{location}'
data4 = f'/media/wencan/Expansion/temp/231124/postprocess/statistics/profile_{location}'
data5 = f'/media/wencan/Expansion/temp/241018/postprocess/statistics/profile_{location}'

datalist = [data0,             data2,             data1,             data4             ,data3                  ]
dy       = [0.0,               0.312,             0.0,               0.312             ,0.07921                ]
color    = ['gray',            'orangered',       'black',           'steelblue'       ,'yellowgreen'          ]
label    = [r'$\mathcal{LS}$', r'$\mathcal{LR}$', r'$\mathcal{HS}$', r'$\mathcal{HR}1$',r'$\mathcal{HR}2$'     ]
lstyle   = ['-',               '-.',              '--',              ':'               ,(0, (3, 1, 1, 1, 1, 1))]
width    = [4.0,               4.0,               4.0,               4.0               ,4.0                    ]
lines = []

# ----------------------------------------------------------------------
# >>> Initialize data                                            ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

for i, datapath in enumerate(datalist):
    
    os.chdir(datapath)
    
    line = ProfileData('profile_mean.dat') 
    
    line.label = label[i]
    line.color = color[i]
    line.width = width[i]
    line.lstyle = lstyle[i]
    
    lines.append(line)

os.chdir( create_folder(outPath) )


# ----------------------------------------------------------------------
# >>> Plot u profile                                             ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_u :
    
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['u'],
                 line.df['y'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.5)
    
    # ax.set_xlim( [-100,507] )
    ax.set_ylim( [0,12] )

    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'u'

    ax.set_ylabel( "$y$", labelpad=-5 )  
    ax.set_xlabel( r'$u$' )
    ax.tick_params( axis='x', pad = 15 )
    ax.tick_params( axis='y', pad = 10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot Reynolds Stress Profile                              ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_RS_uu:

    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )

    for line in lines:
        
        ax.plot( line.df['u`u`'],
                 line.df['y'], 
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    # ax.set_xlim( [0,6000] )
    ax.set_ylim( [0,12] )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)


    # set spacing between major tickers.
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(5.0))

#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_uu'

    ax.set_ylabel( "$y$",labelpad=-5 )  
    ax.set_xlabel( r"$\langle u^{'} u^{'} \rangle $" )
    ax.tick_params( axis='x', pad=15 )
    ax.tick_params( axis='y',pad = 10)
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


if plt_RS_vv:
    
    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )
    
    for line in lines:
        
        ax.plot( line.df['v`v`'], 
                 line.df['y'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
    
    # ax.set_xlim( [0,6000] )
    ax.set_ylim( [0,12] )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5.0))


#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_vv'
    
    ax.set_ylabel( "$y$", labelpad=-5 )  
    ax.set_xlabel( r"$\langle v^{'} v^{'} \rangle$" )
    ax.tick_params( axis='x', pad = 15 )
    ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()
    

if plt_RS_ww:
    
    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )
    
    for line in lines:
        
        ax.plot( line.df['w`w`'], 
                 line.df['y'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
    
    # ax.set_xlim( [0,6000] )
    ax.set_ylim( [0,12] )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)

#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_ww'
    
    ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
    ax.set_ylabel( r"$\rho \langle w^{'} w^{'} \rangle / \tau_w$")
    ax.tick_params( axis='x', pad = 15 )
    ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> plot Reynolds stress                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_RS_uv:

    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )
    
    for line in lines:
        
        ax.plot( -line.df['u`v`'], 
                 line.df['y'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)
    
    # ax.set_xlim( [0,5000] )
    ax.set_ylim( [0,12] )

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)

#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'RS_uv'
    
    ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
    ax.set_ylabel( r"$\rho \langle u^{'} v^{'} \rangle / \tau_w$" )
    ax.tick_params( axis='x', pad = 15 )
    ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> Plot temperature profile                                  ( 5 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_T :
    
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)

    for line in lines:
        
        ax.plot( line.df['T'], 
                 line.df['y'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.5)
    
    # ax.set_xlim( [250,290] )
    ax.set_ylim( [0,12] )
    
    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'T'
    ax.set_ylabel( "$y$", labelpad=-5 )  
    ax.set_xlabel( r'$T$' )
    ax.tick_params( axis='x', pad = 15 )
    ax.tick_params( axis='y', pad = 10 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()

    

# ----------------------------------------------------------------------
# >>> plot total temperature                                     (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/04/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_Tt:
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['Tt'], 
                 line.df['y'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.5)
    
    # ax.set_xlim( [250,290] )
    ax.set_ylim( [0,12] )
    
    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'Tt'
    ax.set_ylabel( "$y$", labelpad=-5 )  
    ax.set_xlabel( r'$T_{t}$' )
    ax.tick_params( axis='x', pad = 15 )
    ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()

# ----------------------------------------------------------------------
# >>> plot total pressure                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/04/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_pt:
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['pt'], 
                 line.df['y'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=10,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=5,
                   width=1.5)
    

    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'pt'
    ax.set_ylabel( "$y$", labelpad=-5 )  
    ax.set_xlabel( r'$p_{t}$' )
    ax.tick_params( axis='x', pad = 15 )
    ax.tick_params( axis='y', pad = 10 )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()


# ----------------------------------------------------------------------
# >>> plot tke                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/04/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_tke:

    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )

    for line in lines:
        
        ax.plot( line.df['u`u`'] + line.df['v`v`'] + line.df['w`w`'],
                 line.df['y'], 
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)

#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.2)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'tke'
    
    ax.set_ylim( [0,12] )
    ax.set_xlabel( "$y_s^+$",labelpad=-5 )  
    ax.set_ylabel( r"$\rho \langle tke \rangle / \tau_w$" )
    ax.tick_params( axis='x', pad=15 )
    ax.tick_params( axis='y',pad = 10)
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()

# ----------------------------------------------------------------------
# >>> plot p`                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/05/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_p_fluc:

    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )

    for line in lines:
        
        ax.plot( line.df['p`']/45447.289, 
                 line.df['y'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)
    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'p_fluc_not_normalized'
    
    ax.set_ylabel( "$y$",labelpad=-5 )  
    ax.set_xlabel( r"$\sqrt{\langle p'p' \rangle} / p_{\infty}$" )
    ax.tick_params( axis='x', pad=15 )
    ax.tick_params( axis='y', pad=10 )
    ax.set_ylim( [0,12] )
#        ax.legend( ) 

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    plt.savefig( figname + fmt )
    plt.show()