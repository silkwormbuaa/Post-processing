#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_profile_grid_study.py
@Time    :   2025/10/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys

import matplotlib
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker
import matplotlib.markers as     markers
import numpy              as     np

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line         import ProfileData
from   vista.directories  import create_folder

plt.rcParams["text.usetex"]         = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family']         = 'Times New Roman'
plt.rcParams['font.size']           = 40

# =============================================================================

plt_u       = False
plt_u_vd    = True
plt_RS      = True

fmt = '.png'

# =============================================================================

OutPath  = '/home/wencan/temp/DataPost/midRe/profile_upstream_grid_study/'

data0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/profile_upstream'
data1 = '/media/wencan/Expansion/temp/smooth_mid/postprocess/statistics/profile_upstream'
data2 = '/media/wencan/Expansion/temp/220927/postprocess/statistics/profile_upstream'
data3 = '/media/wencan/Expansion/temp/241030/postprocess/statistics/profile_upstream'
data4 = '/media/wencan/Expansion/temp/231124/postprocess/statistics/profile_upstream'
data5 = '/media/wencan/Expansion/temp/smooth_mid_x2/postprocess/statistics/profile_upstream'

#data0 = '/home/wencanwu/my_simulation/temp/smooth_wall/x_-53.6.dat'

data250  = source_dir + '/database/Pirozzoli/M2_Retau_250'
data1000 = source_dir + '/database/Pirozzoli/M2_Retau_1000'

datalist = [data0,          data1,         data2,   data3,   data4, data5]
dy       = [0.0,            0.0,           0.312,   0.07921, 0.312, 0.0 ]
color    = ['gray',         'black',       'orangered',  'yellowgreen', 'steelblue', 'black']
label    = ['lowRe_smooth', 'midRe_smooth','lowRe_rough','midRe_0.026', 'midRe_0.1', 'mid_x2']
lstyle   = ['-',            '--',             '-.',      (0, (3, 1, 1, 1, 1, 1)),  ':', ':']
width    = [4.0,            4.0,           4.0,        4.0,        4.0,    4.0]
selected = [ data5 ]
lines    = []

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

for datapath in selected:
    
    i = datalist.index( datapath )
    
    os.chdir(datapath)
    
    line = ProfileData('profile_mean.dat') 
    line.shift_y( dy[i] )
    line.inner_normalize('wall_statistics.dat')
    line.vd_transform()
    
    line.label  = label[i]
    line.color  = color[i]
    line.width  = width[i]
    line.lstyle = lstyle[i]
    
    line.df.to_string( "profile_normalized.dat",
                       index=False, 
                       float_format='%15.7f',
                       justify='left' )
    
    lines.append(line)  

line1000        = ProfileData( data1000 )
line1000.df     = line1000.sparse_log('y+', 0.04)
line1000.label  = "DNS"
line1000.color  = 'black'
line1000.lstyle = 's'
line1000.marker = markers.MarkerStyle(marker='s')

print( OutPath)
os.chdir( create_folder(OutPath) )

def adjust_plot( ax:plt.Axes ):
    
    ax.minorticks_on()
    
    ax.set_xscale( "symlog", linthresh = 1 )

    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=1.5)

    x_minor = ticker.LogLocator( base=10.0, subs = np.arange(1.0,10.0), numticks=100 )
    ax.xaxis.set_minor_locator( x_minor )
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    ax.set_xlabel( "$y_s^+$", labelpad=-5 )  
    ax.tick_params( axis='x', pad=15 )

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)

# ----------------------------------------------------------------------
# >>> plot velocity profile                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/02/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_u:
    
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u+'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.plot( line1000.df['y+'], 
             line1000.df['u+'],
             line1000.color, 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')

    adjust_plot( ax )
    
    ax.set_xlim( [1,3000] )
    ax.set_ylim( [0,25] )
    

    # Adjust the spacing around the plot to remove the white margin
    
    ax.set_ylabel( r'$u^+$' )
    ax.tick_params( axis='y', pad = 10 )

    figname = 'dns_u'

    plt.savefig( figname + fmt )
    print(f"{figname.ljust(15)} is output.")
    plt.close()

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/02/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_u_vd:
    
    fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)

    # law of wall
    
    x1 = np.linspace(1,17,50)
    y1 = np.linspace(1,17,50)
    
    x2 = np.linspace(2,3000,300)
    y2 = np.log(x2)/0.41 + 5.2   

    ax.plot( x1,y1,'gray',ls='-.',lw=4.0)
    ax.plot( x2,y2,'gray',ls='-.',lw=4.0)
    
    text = r"$\frac{\mathrm{log}y^+}{0.41}+5.2$"
    ax.text(2,24,text,fontsize=30)
    ax.plot(3,np.log(3.0)/0.41+5.2,'o',40,color='black')
    ax.plot([3,3],[np.log(3.0)/0.41+5.2,23],'black')
    
    ax.text(13,20,r"$y^+$",fontsize=30)
    ax.plot(14,14,'o',40,color='black')
    ax.plot([14,14],[14,19],'black')    
    
    # profiles
    
    for line in lines:
        
        ax.plot( line.df['ys+'], 
                 line.df['u+_vd'],
                 line.color,   
                 label = line.label, 
                 ls    = line.lstyle,
                 linewidth = line.width)

    ax.plot( line1000.df['y+'], 
             line1000.df['u_vd+'],
             line1000.color, 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')

    adjust_plot( ax )

    ax.set_xlim( [1,3000] )
    ax.set_ylim( [0,28] )
    ax.set_ylabel( r'$\langle u \rangle ^+_{vd}$' )
    ax.tick_params( axis='y', pad=10 )
    
    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'dns_u_vd'

    plt.savefig( figname + fmt )
    print(f"{figname.ljust(15)} is output.")
    plt.close()

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/02/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

if plt_RS:
    
    fig, ax = plt.subplots( figsize=[9,8],constrained_layout=True )

    ax.plot( line1000.df['y+'], 
             line1000.df['urms+']**2*line1000.df['sqrt(rho/rho_w)']**2,
             'black',
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')

    ax.plot( line1000.df['y+'], 
             line1000.df['vrms+']**2*line1000.df['sqrt(rho/rho_w)']**2,
             'black', 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')
    
    ax.plot( line1000.df['y+'], 
             (line1000.df['wrms+'])**2*line1000.df['sqrt(rho/rho_w)']**2,
             'black', 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')
    
    ax.plot( line1000.df['y+'], 
             (line1000.df['uv+'])*line1000.df['sqrt(rho/rho_w)']**2,
             'black', 
             marker='s',
             markersize=10,
             fillstyle='none',
             linestyle='None')
    
    ax.plot( lines[0].df['y+'], 
             lines[0].df['u`u`+'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)

    ax.plot( lines[0].df['y+'], 
             lines[0].df['v`v`+'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)

    ax.plot( lines[0].df['y+'], 
             lines[0].df['w`w`+'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)

    ax.plot( lines[0].df['y+'], 
             lines[0].df['u`v`+'],
             lines[0].color, 
             label = lines[0].label, 
             ls    = lines[0].lstyle,
             linewidth = lines[0].width)
    
    adjust_plot( ax )
    
    ax.set_xlim( [1,3000] )
    ax.set_ylim( [-1.5,9.0] )
    ax.set_ylabel( r"$\langle \rho \rangle \langle u_i'u_j'\rangle /\tau_w$" )
    ax.tick_params( axis='y', pad=10 )
    
    # notations
    
    ax.text( 45,   7, r"$u'u'$",fontsize=30)
    ax.text( 20, 2.5, r"$w'w'$",fontsize=30)
    ax.text( 50, 0.6, r"$v'v'$",fontsize=30)
    ax.text(100,-0.7, r"$u'v'$",fontsize=30)

    # Adjust the spacing around the plot to remove the white margin
    
    figname = 'dns_RS'

    plt.savefig( figname + fmt )
    print(f"{figname.ljust(15)} is output.")
    plt.show()
