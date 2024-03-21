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
width    = [2.0,   2.0,      2.0,    2.0]

plt_pt   = True
plt_Tt   = True
plt_t    = True
plt_tke  = True

# shift y? 
shift_y = True

fmt  = '.png' # or '.pdf'
zoom = True


# ----------------------------------------------------------------------
# >>> define style                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/21  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plt_style( ):
    
    ax = plt.gca()
    
    ax.minorticks_on()
    
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=2.0,
                   pad=10)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=1.5,
                   pad=10)

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)

# =============================================================================

# different locations for all cases
for i in range(6):
    
    lines = []  
    
    # different lines of different cases on one figure
    for j in range(4):
        filename = datapath[j] + f'/profile_mean_{i:02d}.dat'
        line = ProfileData(filename)
        line.shift_y( dy[j] )
        
        line.label = r'$H/\delta_0=$' + label[j]
        line.color = color[j]
        line.width = width[j]
        line.lstyle = lstyle[j]
        
        lines.append(line)
        
    lines[0].label = 'smooth'
    
    os.chdir(outputpath)



# plot total pressure profile
# =============================================================================
    if plt_pt:

        fig, ax = plt.subplots(figsize=[8,8], constrained_layout=True)

        for line in lines:

            x = line.df['pt']/(45447.289*7.824)
            if shift_y: y = line.df['ys']/5.2
            else: y = line.df['y']/5.2
            
            ax.plot(x,y, 
                    label=line.label, 
                    color=line.color, 
                    linestyle=line.lstyle, 
                    linewidth=line.width)
            
            
        ax.set_xlabel(r'$p_t/p_{t_0}$')
        ax.set_xlim(0.2, 1.1)
        
        plt_style()
        
        if shift_y: ax.set_ylabel(r'$y_s/\delta_0$')
        else:       ax.set_ylabel(r'$y/\delta_0$')
        
        if zoom: ax.set_ylim(-0.1,2.0); filename = f'pt_{i}_zoom'
        else:    ax.set_ylim(-0.1,8.0); filename = f'pt_{i}'
        
        if shift_y: filename = filename + '_ys' + fmt
        else:       filename = filename + '_y' + fmt
     
        plt.savefig( filename )
        plt.close()
        

# plot total temperature profile
# =============================================================================

    if plt_Tt:
        
        fig, ax = plt.subplots(figsize=[8,8], constrained_layout=True)

        for line in lines:

            x = line.df['Tt']/288.15                # 160.15
            if shift_y: y = line.df['ys']/5.2
            else: y = line.df['y']/5.2
            
            ax.plot(x,y, 
                    label=line.label, 
                    color=line.color, 
                    linestyle=line.lstyle, 
                    linewidth=line.width)
            

        ax.set_xlabel(r'$T_t/T_{t0}$')
        ax.set_xlim(0.95,1.02)

        plt_style()
         
        if shift_y: ax.set_ylabel(r'$y_s/\delta_0$')
        else:       ax.set_ylabel(r'$y/\delta_0$')
        
        if zoom: ax.set_ylim(-0.1,2.0); filename = f'Tt_{i}_zoom'
        else:    ax.set_ylim(-0.1,8.0); filename = f'Tt_{i}'
        
        if shift_y: filename = filename + '_ys' + fmt
        else:       filename = filename + '_y' + fmt
        
        plt.savefig( filename )
        plt.close()
        
        
# plot temperature profile
# =============================================================================    

    if plt_t:
        
        fig, ax = plt.subplots(figsize=[8,8], constrained_layout=True)

        for line in lines:
        
            x = line.df['T']/160.15
            if shift_y: y = line.df['ys']/5.2
            else: y = line.df['y']/5.2
            
            ax.plot(x,y,
                    label=line.label, 
                    color=line.color, 
                    linestyle=line.lstyle, 
                    linewidth=line.width)
            
            
        ax.set_xlabel(r'$T/T_{\infty}$')
        ax.set_xlim(1.0,1.8)
        
        plt_style()
        
        if shift_y: ax.set_ylabel(r'$y_s/\delta_0$')
        else:       ax.set_ylabel(r'$y/\delta_0$')
        
        if zoom: ax.set_ylim(-0.1,2.0); filename = f'T_{i}_zoom'
        else:    ax.set_ylim(-0.1,8.0); filename = f'T_{i}'
        
        if shift_y: filename = filename + '_ys' + fmt
        else:       filename = filename + '_y' + fmt
        
        plt.savefig( filename )
        plt.close()
        

# plot tke profile
# =============================================================================    

    if plt_tke:
            
        fig, ax = plt.subplots(figsize=[8,8], constrained_layout=True)

        for line in lines:

            x = line.df['tke']
            if shift_y: y = line.df['ys']/5.2
            else: y = line.df['y']/5.2
            
            ax.plot(x,y, 
                    label=line.label, 
                    color=line.color, 
                    linestyle=line.lstyle, 
                    linewidth=line.width)
        
            
        ax.set_xlabel('tke')
        ax.set_xlim(0.0,8000)
        
        plt_style()
        
        if shift_y: ax.set_ylabel(r'$y_s/\delta_0$')
        else:       ax.set_ylabel(r'$y/\delta_0$')
        
        if zoom: ax.set_ylim(-0.1,2.0); filename = f'tke_{i}_zoom'
        else:    ax.set_ylim(-0.1,8.0); filename = f'tke_{i}'
        
        if shift_y: filename = filename + '_ys' + fmt
        else:       filename = filename + '_y' + fmt
        
        plt.savefig( filename )
        plt.close()

    