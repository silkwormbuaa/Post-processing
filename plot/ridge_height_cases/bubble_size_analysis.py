#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size_analysis.py
@Time    :   2024/04/29 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import LineData
from   vista.psd         import pre_multi_psd

file0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/bubble/bubble_size.dat'
file1 = '/media/wencan/Expansion/temp/240211/postprocess/bubble/bubble_size.dat'
file2 = '/media/wencan/Expansion/temp/220927/postprocess/bubble/bubble_size.dat'
file3 = '/media/wencan/Expansion/temp/240210/postprocess/bubble/bubble_size.dat'

outputpath = '/media/wencan/Expansion/temp/DataPost/lowRe_ridge_height/bubble_size'

files = [file0, file1, file2, file3]
tags  = ['smooth_adiabatic', 'H/delta=0.05', 'H/delta=0.1', 'H/delta=0.2'] 
color    = ['gray','black', 'red', 'blue']
lstyle   = ['--',  ':',     '-.',    (0, (3, 1, 1, 1, 1, 1))]

lines = [LineData(file) for file in files]

os.chdir( outputpath )
fig, ax = plt.subplots( figsize=(16,5) )

with open("bubble_size_analysis.dat",'w') as f:
    
    f.write("tag         average_size         st_dev         st_dev_norm\n")

    for i, line in enumerate(lines):
        
        line.df['time'] = np.linspace(20,61,4101)
        
        ax.plot( line.df['time'], line.df['bubble_volume'], color=color[i], linewidth=1.5)
        
        mean_size = np.mean( line.df['bubble_volume'] )
        std_dev = np.std( line.df['bubble_volume'] )
        std_dev_norm = np.std( line.df['bubble_volume'] )/mean_size
        
        print(f"Mean bubble size: {mean_size}, std_dev: {std_dev}, std_dev_norm: {std_dev_norm}")
            
        f.write(f"{tags[i]}         {mean_size}         {std_dev}         {std_dev_norm}\n")

plt.savefig('bubble_size_fluctuation.png')
plt.show()
plt.close()

# pre-multiplied psd of bubble size

fig, ax = plt.subplots()
for i,line in enumerate(lines):

    freq, pm_psd = pre_multi_psd( line.df['bubble_volume'], 100, 10, 0.5 )
    
    ax.semilogx( freq*5.2/507*9.6287, pm_psd, color=color[i], linestyle=lstyle[i], linewidth=1.5 )
    
    
ax.set_xlabel(r'$St_{L_{sep}}$')
ax.set_ylabel(r'$pre-multiplied psd$')
ax.set_ylim( [0.0,0.8] )

plt.savefig(f'psd_bubble.png')
plt.show()
plt.close()