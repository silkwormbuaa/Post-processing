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

from   vista.directories import Directories
from   vista.tools       import find_indices

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

output_dir = '/home/wencan/temp/DataPost/midRe/xy_plane/'
casename   = ['smooth_adiabatic','220927','smooth_mid','231124','241030']
label      = [r'$\mathcal{LS}$',  r'$\mathcal{LR}$', r'$\mathcal{HS}$', r'$\mathcal{HR1}$',r'$\mathcal{HR2}$']
color      = ['gray',             'orangered',       'black',           'steelblue',        'yellowgreen' ]
lstyle     = ['-',                '-.',               '--',             ':',                (0, (3, 1, 1, 1, 1, 1))]

def main():
    
    compute_slope_local()
    
    print("========= spanwise average =========")
    
    compute_slope_rev_spanave()
    
    print("========= dividing streamline =========")
    
    compute_slope_divid_spanave()
    print("========= done =========")


def compute_slope(i,line):
    
    x = np.array(line[:,0])
    y = np.array(line[:,1])
    
    index_max      = np.argmax(y)-10
    index_start, _ = find_indices( y[:index_max], 0.15, mode='sequential')
    
    slope = (y[index_max] - y[index_start]) / (x[index_max] - x[index_start])
    print(f"{casename[i].ljust(20)} slope: {slope:10.4f}, angle: {np.degrees(np.arctan(slope)):10.4f}")
    
    return index_start, index_max
    
def compute_slope_local():
        
    df_dirs  = [f'/media/wencan/Expansion/temp/{case}/postprocess/statistics/xy_planes/' for case in casename]


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
                
                compute_slope(j,line)

        ax.set_xlim([-18,10])
        ax.set_ylim([-0.2, 1])
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(2)
        
        plt.savefig( output_dir + figname + format )
        plt.close()


def compute_slope_rev_spanave():
    
    fig, ax = plt.subplots( figsize=(9,6) )
    
    for i, case in enumerate(casename):
        
        dirs = Directories(f"/home/wencan/temp/{case}")
        os.chdir( dirs.pp_z_average )
        
        with open('sepline.pkl','rb') as f:
            lines = pickle.load(f)
            
            for line in lines:
                
                i_s,i_e = compute_slope(i,line)
                
                ax.plot( line[:,0], line[:,1], color=color[i], 
                         linestyle=lstyle[i], linewidth=1.5 )

                ax.plot( line[i_s:i_e,0], line[i_s:i_e,1], color=color[i], 
                         markersize=3, marker='o', markerfacecolor=color[i] )
                

    ax.set_xlim([-18,10])
    ax.set_ylim([-0.2, 1])
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)
        
    plt.savefig( output_dir + 'spanave_bubble.png')
    plt.close() 


def compute_slope_divid_spanave():
    
    fig, ax = plt.subplots( figsize=(9,6) )
    
    for i, case in enumerate(casename):
        
        dirs = Directories(f"/home/wencan/temp/{case}")
        os.chdir( dirs.pp_z_average )
        
        with open('dividing_streamline.pkl','rb') as f:
            lines = pickle.load(f)
            
            for line in lines:
                
                i_s,i_e = compute_slope(i,line)
                
                ax.plot( line[:,0], line[:,1], color=color[i], 
                         linestyle=lstyle[i], linewidth=1.5 )

                ax.plot( line[i_s:i_e,0], line[i_s:i_e,1], color=color[i], 
                         markersize=3, marker='o', markerfacecolor=color[i] )
                

    ax.set_xlim([-18,10])
    ax.set_ylim([-0.2, 1.5])
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)
        
    plt.savefig( output_dir + 'spanave_divid.png')
    plt.close() 

if __name__ == "__main__":
    main()