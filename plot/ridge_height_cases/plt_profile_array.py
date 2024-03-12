#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_compareprofile.py
@Time    :   2024/03/12 
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

data0 = '/media/wencanwu/Seagate Expansion Drive1/temp/smooth/smooth'
data1 = '/media/wencanwu/Seagate Expansion Drive1/temp/240211/postprocess/statistics/profile_array'
data2 = '/media/wencanwu/Seagate Expansion Drive1/temp/220927/postprocess/statistics/profile_array'
data3 = '/media/wencanwu/Seagate Expansion Drive1/temp/240210/postprocess/statistics/profile_array'
datapath = [data0, data1, data2, data3]

outputpath = '/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe_ridge_height/profile_array'

cases = ['smooth', '240211', '220927', '240210']

dy       = [0,     0.156,   0.312,        0.624]
color    = ['gray','black', 'red', 'blue']
label    = ['',    '0.05',  '0.10','0.20']
lstyle   = ['--',  ':',     '-.',    (0, (3, 1, 1, 1, 1, 1))]
width    = [4.0,   4.0,      4.0,    4.0]

# different locations for all cases
for i in range(6):
    
    lines = []  
    
    # different lines of different cases on one figure
    for j in range(4):
        filename = datapath[j] + f'/profile_mean_{i:02d}.dat'
        line = ProfileData(filename)
        line.shift_y( dy[j] )
        line.df['y_'] = line.df['y']/5.2
        
        line.label = r'$H/\delta_0=$' + label[j]
        line.color = color[j]
        line.width = width[j]
        line.lstyle = lstyle[j]
        
        lines.append(line)
        
    lines[0].label = 'smooth'
    
    os.chdir(outputpath)

    fig, ax = plt.subplots(figsize=[12,15], constrained_layout=True)

    for line in lines:
        
        ax.plot(line.df['pt']/45447.289, 
                line.df['ys']/5.2,  
                label=line.label, 
                color=line.color, 
                linestyle=line.lstyle, 
                linewidth=line.width)
        
        ax.minorticks_on()
    
    ax.set_ylim(-0.1,8.0)
#    ax.set_xlim(0, 1.2)
    fig.legend(loc='upper left')
    plt.savefig(f'ys_pt_{i}.png')
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
        
    ax.set_ylim(-0.1,8.0)
#    ax.set_xlim(0, 1.2)
    fig.legend(loc='upper left')
    plt.savefig(f'ys_Tt_{i}.png')
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
    
    ax.set_ylim(-0.1,8.0)
#    ax.set_xlim(0, 1.2)
    fig.legend(loc='upper left')
    plt.savefig(f'ys_T_{i}.png')
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
    
    ax.set_ylim(-0.1,8.0)
#    ax.set_xlim(0, 1.2)
    fig.legend(loc='upper left')
    plt.savefig(f'ys_tke_{i}.png')
#    plt.show()
    plt.close()
    