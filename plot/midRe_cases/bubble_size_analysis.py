#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size_analysis.py
@Time    :   2024/09/09
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plot the bubble size fluctuation and pre-multiplied psd
'''

import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt
from   scipy.stats       import skew, kurtosis

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import LineData
from   vista.psd         import pre_multi_psd
from   vista.directories import create_folder

# =============================================================================

file       = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/bubble/bubble_size.dat'
outputpath = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/bubble'
lsep       = 49.504  # 47.944
delta0     = 5.2

files = [file]
tags  = ['lowRe_smooth'] 
color    = ['blue']
lstyle   = ['-']

# =============================================================================

lines = [LineData(file) for file in files]

os.chdir( create_folder(outputpath) )
fig, ax = plt.subplots( figsize=(16,5) )

with open("bubble_size_analysis.dat",'w') as f:
    
    f.write("tag         average_size         st_dev         st_dev_norm\n")

    for i, line in enumerate(lines):
        
        line.df['time'] = np.linspace(20,61,4101)
        mean_size = np.mean( line.df['bubble_volume'] )
        std_dev = np.std( line.df['bubble_volume'] )
        std_dev_norm = np.std( line.df['bubble_volume'] )/mean_size
        
        ax.plot( (np.array(line.df['itime'])- 20.0 ) * 507.0 * (1.0/delta0), 
                 (np.array(line.df['bubble_volume'])-mean_size)/mean_size, 
                 color=color[i], 
                 linewidth=1.5)
        
        
        print(f"Mean bubble size: {mean_size}, std_dev: {std_dev}, std_dev_norm: {std_dev_norm}")
            
        f.write(f"{tags[i]}         {mean_size}         {std_dev}         {std_dev_norm}\n")

plt.savefig('bubble_size_fluctuation.png')
plt.show()
plt.close()

# pre-multiplied psd of bubble size

fig, ax = plt.subplots(figsize=(12,6))
for i,line in enumerate(lines):

    freq, pm_psd = pre_multi_psd( (np.array(line.df['bubble_volume'])-mean_size)/mean_size, 
                                  100, 8, 0.5, 
                                  nfft=len(line.df['bubble_volume']) )
    
    ax.semilogx( freq*lsep/507, pm_psd, color=color[i], linestyle=lstyle[i], linewidth=1.5 )
    
    
ax.set_xlabel(r'$St_{L_{sep}}$')
ax.set_ylabel(r'$pre-multiplied psd$')
ax.set_ylim( [0.0,0.8] )

plt.savefig(f'psd_bubble.png')
plt.show()
plt.close()

# =============================================================================
# skewness and flatness

for line in lines:
    
    mean_size = np.mean( line.df['bubble_volume'] )
    bubble_vol_fluc = (np.array(line.df['bubble_volume'])-mean_size)/mean_size
    skewness = skew( bubble_vol_fluc )
    flatness = kurtosis( bubble_vol_fluc )
    
    print(f"Skewness: {skewness}, Flatness: {flatness}")