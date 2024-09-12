#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   compare_bubble_shape.py
@Time    :   2024/05/06 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import pickle
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

output_dir = '/media/wencan/Expansion/temp/DataPost/lowRe_ridge_height/xy_plane/'

casename = ['smooth_adiabatic','240211','220927','240210']

df_dirs  = [f'/media/wencan/Expansion/temp/{case}/postprocess/statistics/xy_planes/' for case in casename]

label    = ['smooth',  r'$H/\delta_0=0.05$', r'$H/\delta_0=0.10$', r'$H/\delta_0=0.20$',]
color    = ['gray','black', 'red', 'blue']
lstyle   = ['--',  ':',     '-.',    (0, (3, 1, 1, 1, 1, 1))]

figname0 = 'bubble_shape'
format   = '.png'

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

    ax.set_xlim([-18,10])
    ax.set_ylim([-0.2, 1])
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)
    
    plt.savefig( output_dir + figname + format )
    plt.close()
    
    