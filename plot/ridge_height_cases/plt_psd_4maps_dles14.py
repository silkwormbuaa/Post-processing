#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_psd_4maps.py
@Time    :   2024/05/20
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt
from   matplotlib.colors import Normalize
from   matplotlib.gridspec import GridSpec

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.tools       import get_filelist

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

figname = 'psd_4maps_dles14'
fmt = '.pdf'

outpath = '/media/wencan/Expansion/temp/DataPost/lowRe_ridge_height/psd'

datapaths = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/probes/psd_ridge'
datapath1 = '/media/wencan/Expansion/temp/240211/postprocess/probes/psd_ridge'
datapath2 = '/media/wencan/Expansion/temp/220927/postprocess/probes/psd_ridge'
datapath3 = '/media/wencan/Expansion/temp/240210/postprocess/probes/psd_ridge'
datapaths = [datapaths, datapath1, datapath2, datapath3]

# separation range in normalized unit
x_separs = [-8.42  , -9.510, -10.65744, -13.85 ]
x_reatts = [1.06795,   1.93,   2.468831,  3.25 ]
x_ppmaxs = [-7.3333,  -8.7426, -10.0611,  -11.5126 ]

# labels
label    = ['smooth',  r'$H/\delta_0=0.05$', r'$H/\delta_0=0.10$', r'$H/\delta_0=0.20$',]

# subplots position and size

fig = plt.figure( figsize=(30, 12.5) )
gs = GridSpec(2,20)
ax_ranges = [0,11,0,10]
ax_rangee = [9,20,9,20]
colornorm = Normalize( vmin=0, vmax=0.3 )
axs = list()

for i in range(4):
    
    psdfilelist = get_filelist( datapaths[i] )
    
    x = [] ; st = [] ; pmpsd = []
    
    for j, psdfile in enumerate( psdfilelist ):
        
        print(f"Imported {j+1:05d} probes: {psdfile[-20:]}")
        
        probe = ProbeData( )
        probe.read_psd( psdfile )
        st_probe = np.array( probe.psd_df['freq'] )*5.2/507.0
        x_probe  = [(probe.xyz[0]-50.4)/5.2]*len(st_probe)
        pmpsd_probe = np.array( probe.psd_df['pmpsd_p_fluc'] )

        x.append( x_probe )
        st.append( st_probe )
        pmpsd.append( pmpsd_probe )
    
    st_sep = np.array( st ) * (x_reatts[i]-x_separs[i])
    
    ax = fig.add_subplot( gs[i//2,ax_ranges[i]:ax_rangee[i]] )
    pmpsd = ax.pcolormesh(x,st_sep,pmpsd,
                          shading='gouraud',
                          cmap='Greys',
                          norm=colornorm)
    
    ax.plot([x_separs[i], x_separs[i]], [0.001, 100], 'black', lw=2.5, linestyle='--')
    ax.plot([x_reatts[i], x_reatts[i]], [0.001, 100], 'black', lw=2.5, linestyle='--')
    ax.plot([x_ppmaxs[i], x_ppmaxs[i]], [0.001, 100], 'black', lw=2.5)
    
    ax.text( 4.0, 0.02, label[i] )
    
    ax.set_xlim([-12.5, 10.0])
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
    
    ymajor_ticks = np.array([0.01, 0.1, 1, 10, 100])
    minor_ticks = []
    
    for j in range(len(ymajor_ticks)-1):
        minor_ticks.extend(np.linspace(ymajor_ticks[j], ymajor_ticks[j+1],num=9)[1:-1])

    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(minor_ticks, minor=True)
    
    if i % 2 == 0:
        ylabel = r'$St = fL_{sep}/u_{\infty}$'
        ax.set_ylabel(ylabel)
    
    if i < 2:
        ax.xaxis.set_ticklabels([])
    else:
        xlabel = r'$(x-x_{imp})/\delta_0$'
        ax.set_xlabel(xlabel)
    
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)
    
    axs.append(ax)
    
        
axs[3].set_xlim([-15.0, 10.0])

cbar_label = r'$f \cdot \mathrm{PSD}(f)/ \int \mathrm{PSD}(f) \mathrm{d} f$'
cbar_ticks = np.linspace(0, 0.3, 7)
#cbar_ax = fig.add_subplot([0.4,0.05,0.6,0.1],visible=False)
#cbar_ax.set_position([0.0, 0.0, 0.7, 0.03])
#cbar = fig.colorbar(pmpsd, ax=cbar_ax, orientation='horizontal', shrink=0.8,
#                    ticks=cbar_ticks)
#cbar.ax.tick_params(axis='x', direction='in', bottom=True, top=False,
#                    length=10.0, width=1.5)
#cbar.ax.set_ylabel(cbar_label, rotation='horizontal', labelpad=0.2)
#cbar.ax.yaxis.set_label_coords(-0.45,-0.15)


plt.subplots_adjust(left=0.07, right=0.99, bottom=0.15, top=0.95, wspace=0.2, hspace=0.2)
os.chdir( outpath )
plt.savefig( figname + fmt )