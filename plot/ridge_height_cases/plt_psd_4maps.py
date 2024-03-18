#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_psd_4maps.py
@Time    :   2024/03/18 
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

from   vista.probe       import ProbeData
from   vista.tools       import get_filelist

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = "Times New Roman"

outpath = '/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe_ridge_height/psd'

datapaths = '/media/wencanwu/Seagate Expansion Drive1/temp/smooth/postprocess/probes/psd_ridge'
datapath1 = '/media/wencanwu/Seagate Expansion Drive1/temp/240211/postprocess/probes/psd_ridge'
datapath2 = '/media/wencanwu/Seagate Expansion Drive1/temp/220927/postprocess/probes/psd_ridge'
datapath3 = '/media/wencanwu/Seagate Expansion Drive1/temp/240210/postprocess/probes/psd_ridge'
datapaths = [datapaths, datapath1, datapath2, datapath3]

# separation range in normalized unit
x_separs = [-8.560779, -9.510, -10.65744, -13.85 ]
x_reatts = [1.06795,   1.93,   2.468831,  3.25 ]
x_ppmaxs = [-7.3333,  -8.7426, -10.0611,  -11.5126 ]


fig, axs = plt.subplots(4, 1, figsize=(10, 15), sharex=True)

for i, ax in enumerate(axs.flat):
    
    psdfilelist = get_filelist( datapaths[i] )
    
    x = [] ; st = [] ; pmpsd = []
    
    for j, psdfile in enumerate( psdfilelist ):
        
        print(f"Imported {j:05d} probes: {psdfile[-20:]}")
        
        probe = ProbeData( )
        probe.read_psd( psdfile )
        st_probe = np.array( probe.psd_df['freq'] )*5.2/507.0
        x_probe  = [(probe.xyz[0]-50.4)/5.2]*len(st_probe)
        pmpsd_probe = np.array( probe.psd_df['pmpsd_p_fluc'] )

        x.append( x_probe )
        st.append( st_probe )
        pmpsd.append( pmpsd_probe )
    
    st_sep = np.array( st ) * (x_reatts[i]-x_separs[i])
    
    pmpsd = ax.pcolormesh(x,st_sep,pmpsd,
                             shading='gouraud',
                             cmap='Greys',
                             vmin=0, vmax=0.3)
    
    ax.plot([x_separs[i], x_separs[i]], [0.001, 100], 'r', lw=1.0)
    ax.plot([x_reatts[i], x_reatts[i]], [0.001, 100], 'r', lw=1.0)
    ax.plot([x_ppmaxs[i], x_ppmaxs[i]], [0.001, 100], 'b', lw=1.0)
    
    ax.set_xlim([-13.2, 12.0])
    ax.set_yscale('log')
    ax.set_ylim([0.01, 100])
    ax.minorticks_on()
    ax.tick_params(which='major',
                    axis='both',
                    direction='out',
                    length=20,
                    width=2.0)
    ax.tick_params(which='minor',
                    axis='both', 
                    direction='out',
                    length=10,
                    width=1.5)
    
        
axs[3].set_xlim([-15.0, 12.0])

cbar_label = r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$'
cbar_ticks = np.linspace(0, 0.3, 7)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
cbar = plt.colorbar(pmpsd, ax=cbar_ax, orientation='horizontal', 
                    ticks=cbar_ticks)
cbar.ax.tick_params(axis='x', direction='in', bottom=True, top=False,
                    length=15.0, width=1.5)
cbar.ax.set_ylabel(cbar_label, rotation='horizontal', labelpad=0.2)
cbar.ax.yaxis.set_label_coords(-0.4,-0.15)


plt.subplots_adjust( bottom=0.3,right=0.80)

os.chdir( outpath )
plt.savefig( 'psd_4maps.png' )