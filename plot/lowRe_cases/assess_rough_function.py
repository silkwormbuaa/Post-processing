#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   assess_rough_function.py
@Time    :   2024/06/20 
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
import pandas             as     pd

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line         import ProfileData
from   vista.directories  import create_folder
from   vista.tools        import create_linear_interpolator

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 40

fmt = '.png'

OutPath  = '/home/wencanwu/my_simulation/temp/DataPost/profile/rf/'

data0 = '/media/wencanwu/Seagate Expansion Drive1/temp/smooth_isothermal/results/profile'
data1 = '/media/wencanwu/Seagate Expansion Drive1/temp/221014/postprocess/statistics/upstream_profile'
data2 = '/media/wencanwu/Seagate Expansion Drive1/temp/220926/postprocess/statistics/upstream_profile'
data3 = '/media/wencanwu/Seagate Expansion Drive1/temp/220825/postprocess/statistics/upstream_profile'
data4 = '/media/wencanwu/Seagate Expansion Drive1/temp/220927/postprocess/statistics/upstream_profile'
data5 = '/media/wencanwu/Seagate Expansion Drive1/temp/221221/postprocess/statistics/upstream_profile'
data6 = '/media/wencanwu/Seagate Expansion Drive1/temp/smooth_adiabatic/postprocess/statistics/upstream_profile'

# option zone
# =============================================================================

# opt = 1, shift wall to valley
# opt = 2, shift wall coordinate to roughness elevation height
# opt = 3, do not shift, origin at the top

opt   = 3
scale = 1    # 1 inner scale; 2 outer scale
var   = 'u+_vd'  # 'u+' 'u+_vd'

if   var == 'u':
    ylim = [0,1.1]
    ylabel = r'$\langle u \rangle /u_{\infty}$'
else:
    ylim = [0,23]
    if var == 'u+': ylabel = r'$\langle u \rangle ^+$'
    elif var == 'u+_vd': ylabel = r'$\langle u \rangle ^+_{vD}$'


if   scale == 1: 
    figname = var + '_in'
    xlim = [1,1200]; xlim_rf = [20,400]
    xlabel = "$y_s^+$"
    
elif scale == 2: 
    figname = var + '_out'
    xlim = [0,1.5]; xlim_rf = [0,1.2]
    xlabel = "$y_s/\delta_0$"
    
if   opt == 1: 
    dy = [0.0, 0.52,  0.52,  0.52,  0.52,  0.52, 0.0]
    figname += '_valley'
elif opt == 2: 
    dy = [0.0, 0.494, 0.468, 0.416, 0.312, 0.26, 0.0]
    figname += '_elevation'
elif opt == 3: 
    dy = [0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0]
    figname += '_top'


datalist = [data6, data1,   data2,   data3,                   data4,        data5]
color    = ['gray','green', 'blue', 'black',                  'red',        'purple', 'yellow']
label    = ['',    '2.0',   '1.0',  '0.5',                    '0.25',       '0.125', 'smooth']
lstyle   = ['--',  ':',     '-.',    (0, (3, 1, 1, 1, 1, 1)), (0, (10, 3)), '-',     'dotted']
width    = [4.0,   4.0,      4.0,    4.0,                     4.0,          4.0,     4.0]

# =============================================================================


lines = []

for i, datapath in enumerate(datalist):
    
    os.chdir(datapath)
    
    line = ProfileData('profile_mean.dat') 
    line.shift_y( dy[i] )
    line.inner_normalize('wall_statistics.dat')
    line.vd_transform()
    
    line.label = r'$D/\delta_0=$' + label[i]
    line.color = color[i]
    line.width = width[i]
    line.lstyle = lstyle[i]
    
    line.df.to_string( "profile_normalized.dat",
                       index=False, 
                       float_format='%15.7f',
                       justify='left' )
    
    lines.append(line)
    
lines[0].label = 'smooth'

create_folder( OutPath )
os.chdir( OutPath )

fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)

for line in lines:
    
    if scale == 1: 
        wallnorm = line.df['ys+']
        varnorm  = line.df[var]
        
    elif scale == 2: 
        wallnorm = line.df['ys']/5.2
        varnorm  = line.df[var]/507.0
    
    ax.plot( wallnorm, 
             varnorm,
             line.color,   
             label = line.label, 
             ls    = line.lstyle,
             linewidth = line.width)

ax.minorticks_on()

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

if scale == 1:   
    ax.set_xscale( "symlog", linthresh = 1 )
#    x_minor = matplotlib.ticker.LogLocator( base=10.0, subs = np.arange(0.1,10.0) )
#    ax.xaxis.set_minor_locator( x_minor )

# Adjust the spacing around the plot to remove the white margin

ax.set_xlim( xlim )
ax.set_ylim( ylim )
ax.set_xlabel( xlabel, labelpad=-5 )  
ax.set_ylabel( ylabel )
ax.tick_params( axis='x', pad=15 )
ax.tick_params( axis='y', pad=10 )
#        ax.legend( ) 

# set the bounding box of axes
ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(3)

plt.savefig( figname + fmt )
plt.show()


# exit the program if using u 

if var == 'u': sys.exit()

# ----------------------------------------------------------------------
# >>> compute roughness function                                 (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/06/20  - created
#
# Desc
#
# ----------------------------------------------------------------------

ys = np.linspace(400,1200,381)
uvd_int = []
rfs = []

for line in lines:
    
    yst  = line.df['ys+']
    uvdt = line.df['u+_vd']
    
    interpolator = create_linear_interpolator( yst, uvdt )
    
    uvd_int.append( np.array([interpolator(y) for y in ys]) )
    
for i in range(1,len(lines)):
    
    rfs.append( (uvd_int[0] - uvd_int[i]) )
    

# compute the average value and error bar

rf_avg = np.mean( rfs, axis=1 )
rf_err = [np.max(rf)-np.min(rf) for rf in rfs]

print( "uvd average value in the wake: smooth, 2.0, 1.0, 0.5, 0.25, 0.125")
print( np.mean( uvd_int, axis=1) )

print( "\n Roughness function in the wake: 2.0, 1.0, 0.5, 0.25, 0.125")
print( rf_avg )

print( "\n Roughness function error bar in the wake: 2.0, 1.0, 0.5, 0.25, 0.125")
print( rf_err )


# ----------------------------------------------------------------------
# >>> plot roughness function                                    (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/06/20  - created
#
# Desc
#
# ----------------------------------------------------------------------

    
fig, ax = plt.subplots(figsize=[9,8],constrained_layout=True)

for i,rf in enumerate(rfs):
    
    ax.plot( ys, rf,
             lines[i+1].color,
             label = lines[i+1].label,
             ls    = lines[i+1].lstyle,
             linewidth = lines[i+1].width )
    
ax.minorticks_on()
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

if scale == 1:
    ax.set_xlim( [400,1200] )
    ax.set_ylim( [-2,4] )
    xlabel  = "$y_s^+$"

elif scale == 2:
    ax.set_xlim( [0,1.5] )
    ax.set_ylim( [-2,4] )
    xlabel = "$y_s/\delta_0$"
    
# ax.set_xscale( "symlog", linthresh = 1 )
# x_minor = matplotlib.ticker.LogLocator( 
#                     base=10.0, subs = np.arange(1.0,10.0) )
# ax.xaxis.set_minor_locator( x_minor )

# Adjust the spacing around the plot to remove the white margin

figname += '_rf_wake'

ax.set_xlabel( xlabel, labelpad=-5 )  
ax.set_ylabel( r'$\Delta \langle u \rangle ^+_{vD}$' )
ax.tick_params( axis='x', pad=15 )
ax.tick_params( axis='y', pad=10 )

# set the bounding box of axes
ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(3)
plt.savefig( figname + fmt )
plt.show()
