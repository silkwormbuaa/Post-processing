#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   compare_bubble_shape.py
@Time    :   2025/03/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import pickle
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.tools       import find_indices

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

output_dir = '/home/wencan/temp/DataPost/midRe/xy_plane/'

casename = ['smooth_adiabatic','220927','smooth_mid','231124','241030']

df_dirs  = [f'/media/wencan/Expansion/temp/{case}/postprocess/statistics/xy_planes/' for case in casename]

label    = [r'$\mathcal{LS}$',  r'$\mathcal{LR}$', r'$\mathcal{HS}$', r'$\mathcal{HR1}$',r'$\mathcal{HR2}$']
color    = ['gray',             'orangered',       'black',           'steelblue',        'yellowgreen' ]
lstyle   = ['-',                '-.',               '--',             ':',                (0, (3, 1, 1, 1, 1, 1))]

figname0 = 'bubble_shape'
format   = '.png'


def compute_slope(i,line):
    
    x = np.array(line[:,0])
    y = np.array(line[:,1])
    
    index_max      = np.argmax(y)
    index_start, _ = find_indices( y[:index_max], 0.1, mode='sequential')
    
    slope = (y[index_max] - y[index_start]) / (x[index_max] - x[index_start])
    print(f"{casename[i]} slope: {slope}")
    
    
for i in range(2):
    
    if i == 0:   figname = figname0 + '_ridge'
    elif i == 1: figname = figname0 + '_valley'
    
    # create figure
    
    fig, ax = plt.subplots( figsize=(9,6) )
    
    seplinefiles = [df_dir + f'seplines_{i:02d}.pkl' for df_dir in df_dirs]
    
    for j,seplinefile in enumerate(seplinefiles):
        
        with open(seplinefile,'rb') as f: lines = pickle.load(f)
        for line in lines:
            ax.plot( line[:,0], line[:,1], color=color[j], 
                     linestyle=lstyle[j], linewidth=1.5 )
            
            compute_slope(j,line)

    ax.set_xlim([-18,10])
    ax.set_ylim([-0.2, 1])
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)
    
    plt.savefig( output_dir + figname + format )
    plt.close()
    
