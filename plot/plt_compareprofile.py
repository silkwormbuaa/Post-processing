#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_compareprofile.py
@Time    :   2024/02/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   compare profiles at 
'''

import os
import sys
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import ProfileData

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 40

datapath = '/home/wencanwu/my_simulation/temp/DataPost/profile_arrays/'

cases = ['smooth', '221014', '220926', '220825', '220927',  '221221']

dy       = [0,     0.494,   0.468,   0.416,                   0.312,        0.26]
color    = ['gray','green', 'blue', 'black',                  'red',        'purple']
label    = ['',    '2.0',   '1.0',  '0.5',                    '0.25',       '0.125']
lstyle   = ['--',  ':',     '-.',    (0, (3, 1, 1, 1, 1, 1)), (0, (10, 3)), '-']
width    = [4.0,   4.0,      4.0,    4.0,                     4.0,          4.0]

# different locations for all cases
for i in range(1, 7):
    
    lines = []  
    
    # different lines of different cases on one figure
    for j in range(6):
        filename = datapath + cases[j] + f'/profile_mean_{i}.dat'
        line = ProfileData(filename)
        line.shift_y( dy[j] )
        line.df['y_'] = line.df['y']/5.2
        
        line.label = r'$D/\delta_0=$' + label[j]
        line.color = color[j]
        line.width = width[j]
        line.lstyle = lstyle[j]
        
        lines.append(line)
        
    lines[0].label = 'smooth'
    
    os.chdir(datapath)

    fig, ax = plt.subplots(figsize=[12,15], constrained_layout=True)

    for line in lines:
        
        ax.plot(line.df['pt']/45447.289, 
                line.df['ys']/5.2,  
                label=line.label, 
                color=line.color, 
                linestyle=line.lstyle, 
                linewidth=line.width)
        
        ax.minorticks_on()
    
    ax.set_ylim(-0.1,2.0)
#    ax.set_xlim(0, 1.2)
    fig.legend(loc='upper left')
    plt.savefig(f'ys_bottom_pt_{i}.png')
#    plt.show()
    plt.close()

# =============================================================================

    fig, ax = plt.subplots(figsize=[12,15], constrained_layout=True)

    for line in lines:
        
        ax.plot(line.df['Tt']/288.15,    #160.15
                line.df['ys']/5.2,  
                label=line.label, 
                color=line.color, 
                linestyle=line.lstyle, 
                linewidth=line.width)
        
        ax.minorticks_on()
        
    ax.set_ylim(-0.1,2.0)
#    ax.set_xlim(0, 1.2)
    fig.legend(loc='upper left')
    plt.savefig(f'ys_bottom_Tt_{i}.png')
#    plt.show()
    plt.close()
# =============================================================================    

    fig, ax = plt.subplots(figsize=[12,15], constrained_layout=True)

    for line in lines:
        
        ax.plot(line.df['T'], 
                line.df['ys']/5.2,  
                label=line.label, 
                color=line.color, 
                linestyle=line.lstyle, 
                linewidth=line.width)
        
        ax.minorticks_on()
    
    ax.set_ylim(-0.1,2.0)
#    ax.set_xlim(0, 1.2)
    fig.legend(loc='upper left')
    plt.savefig(f'ys_bottom_T_{i}.png')
#    plt.show()
    plt.close()
    
# =============================================================================    

    fig, ax = plt.subplots(figsize=[12,15], constrained_layout=True)

    for line in lines:
        
        ax.plot(line.df['tke'], 
                line.df['ys']/5.2,  
                label=line.label, 
                color=line.color, 
                linestyle=line.lstyle, 
                linewidth=line.width)
        
        ax.minorticks_on()
    
    ax.set_ylim(-0.1,2.0)
#    ax.set_xlim(0, 1.2)
    fig.legend(loc='upper left')
    plt.savefig(f'ys_bottom_tke_{i}.png')
#    plt.show()
    plt.close()
    