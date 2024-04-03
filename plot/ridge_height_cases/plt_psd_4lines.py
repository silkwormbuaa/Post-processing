#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_psd_4lines.py
@Time    :   2024/04/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import pandas            as     pd
import numpy             as     np
from   scipy.interpolate import make_interp_spline

source_dir = os.path.realpath(__file__).split('plot')[0] 
sys.path.append( source_dir )

import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker
from   vista.probe       import ProbeData
from   vista.timer       import timer
from   vista.line        import data_bin_avg

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 30

# option zone
# =============================================================================

loc = 'sep'
independent_len = False
figname = 'psd_4lines'
fmat = '.png'
cases = [0,1,2,3]
withT = [True, True, False, True]
color  = ['gray','black','red','blue']
lstyle = ['--',  ':',     '-.',    (0, (3, 1, 1, 1, 1, 1)) ]
width  = [4.0,   4.0,      4.0,    4.0 ]
label  = ['smooth', r'H/\delta_0=0.05', r'H/\delta_0=0.10', r'H/\delta_0=0.20']
showlegend = False

# =============================================================================

datapath0 = '/media/wencanwu/Seagate Expansion Drive1/temp/smooth_adiabatic/probes'
datapath1 = '/media/wencanwu/Seagate Expansion Drive1/temp/240211/probes/'
datapath2 = '/media/wencanwu/Seagate Expansion Drive1/temp/220927/probes/'
datapath3 = '/media/wencanwu/Seagate Expansion Drive1/temp/240210/probes/'

outpath  = '/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe_ridge_height/psd'

datapaths = [datapath0, datapath1, datapath2, datapath3]

prbfiles_sep   = ['probe_00142.dat', 'probe_00049.dat', 'probe_00118.dat', 'probe_00020.dat']
prbfiles_pfmax = ['probe_00156.dat', 'probe_00058.dat', 'probe_00125.dat', 'probe_00048.dat']  

if loc == 'sep':      prbfiles = prbfiles_sep
elif loc == 'pf_max': prbfiles = prbfiles_pfmax

len_sep = [ 9.6287, 11.44, 13.12627, 17.1 ]

probes = []

with timer("reading probes"):
    
    for casenr in cases:
        
        os.chdir( datapaths[casenr] )
        probe = ProbeData( prbfiles[casenr], withT=withT[casenr] )
        probe.cleandata( t_start=20.0 )
        probe.get_fluc( 'p' )
        probe.pre_multi_psd( 'p_fluc', n_seg=8, overlap=0.5 )
        probes.append( probe )

# =============================================================================
# plot
# =============================================================================

os.chdir( outpath )

fig, ax = plt.subplots( figsize=(9, 8), constrained_layout=True )

for casenr in cases:
    
    probe = probes[casenr]
    
    if independent_len: lsep = len_sep[casenr]
    else:               lsep = len_sep[0]
    
    st = probe.psd_df['freq'] * 5.2 / 507 * lsep
    
    ax.semilogx( st, probe.psd_df['pmpsd_p_fluc'],
                 color=color[casenr], 
                 linestyle=lstyle[casenr], 
                 linewidth=width[casenr], 
                 label=label[casenr] )
    
ax.minorticks_on()
ax.tick_params( which='major',
                axis='both',
                direction='in',
                length=15,
                width=2,
                pad=15)
ax.tick_params( which='minor',
                axis='both', 
                direction='in',
                length=10,
                width=1)

#    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
#    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))

# ax.set_xlim( [0.01,1] )
# ax.set_ylim( [0.0,0.3] )

ax.set_xlabel( r'$f L_{sep}/u_{\infty}$' )
ax.set_ylabel( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$' )

ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(2)

if showlegend: ax.legend()

plt.savefig( figname + fmat )
plt.close()
    
# =============================================================================
# plot bin chart
# =============================================================================

fig, ax = plt.subplots( figsize=(9, 8), constrained_layout=True )

for casenr in cases:
    
    probe = probes[casenr]
    
    if independent_len: lsep = len_sep[casenr]
    else:               lsep = len_sep[0]
    
    st = probe.psd_df['freq'] * 5.2 / 507 * lsep
    
    df = pd.DataFrame( { 'st': st, 'psd': probe.psd_df['pmpsd_p_fluc'] } )
    bin_edges = np.logspace( -2, 2, 41, endpoint=True )
    
    df_bin = data_bin_avg( df, 'st', bin_edges )
    
    # clear data
    
    df_bin.dropna( inplace=True )
    
    # plot    
    
    barwidth = np.diff( bin_edges )
    
    ax.plot( df_bin['st_mid'], df_bin['psd'], 
             color=color[casenr], 
             linewidth=width[casenr], 
             linestyle=lstyle[casenr],
             label=label[casenr] )

ax.set_xscale( 'log' )
ax.set_xlim( [0.01,100] )
ax.set_ylim( [0.0, 0.5] )

ax.minorticks_on()
ax.tick_params( which='major',
                axis='both',
                direction='in',
                length=15,
                width=2,
                pad=15)
ax.tick_params( which='minor',
                axis='both', 
                direction='in',
                length=10,
                width=1)

ax.set_xlabel( r'$f L_{sep}/u_{\infty}$' )
ax.set_ylabel( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$' )

ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(2)

if showlegend: ax.legend()

plt.savefig( figname + '_bin' + fmat )
plt.close()    
    