#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_psd_lines.py
@Time    :   2025/03/04 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import pandas            as     pd
import numpy             as     np

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

loc     = 'psd_max'
independent_len = True
figname = 'psdlines'
fmat    = '.png'
cases   = [0,1,2,3]
withT   = [True, True, True, True]
color   = ['black','yellowgreen','steelblue','firebrick']
lstyle  = ['--', (0, (3, 1, 1, 1, 1, 1)), ':', '-' ]
width   = [4.0,  4.0, 4.0, 4.0 ]
label   = ['low smooth', 'R1', 'R2', 'R3']
showlegend  = True
premultiply = True
normalize   = True

# =============================================================================

datapath0 = '/home/wencan/temp/smooth_adiabatic/probes/'
datapath1 = '/home/wencan/temp/250120/probes/'
datapath2 = '/home/wencan/temp/250218/probes/'
datapath3 = '/home/wencan/temp/250304/probes/'

outpath  = '/home/wencan/temp/DataPost/lowRe/psd'

datapaths = [datapath0, datapath1, datapath2, datapath3]

prbfiles_sep     = ['probe_00144.dat', 'probe_00084.dat', 'probe_00153.dat', 'probe_00133.dat']
prbfiles_pfmax   = ['probe_00158.dat', 'probe_00073.dat', 'probe_00138.dat', 'probe_00145.dat']  

prbfiles_psdmax  = ['probe_00144.dat', 'probe_00073.dat', 'probe_00123.dat', 'probe_00133.dat']

if loc   == 'sep':     prbfiles = prbfiles_sep
elif loc == 'pf_max':  prbfiles = prbfiles_pfmax
elif loc == 'psd_max': prbfiles = prbfiles_psdmax

len_sep = [ 9.52, 16.96, 9.1, 10.56]

probes = []

with timer("reading probes"):
    
    for casenr in cases:
        
        os.chdir( datapaths[casenr] )
        probe = ProbeData( prbfiles[casenr], withT=withT[casenr] )
        probe.cleandata( t_start=20.0 )
        probe.get_fluc( 'p' )
        
        if premultiply: 
            probe.pre_multi_psd( 'p_fluc', n_seg=8, overlap=0.5, normalize=normalize )
            var = 'pmpsd_p_fluc'
        else:         
            probe.compute_psd( 'p_fluc', n_seg=8, overlap=0.5 )
            var = 'psd_p_fluc'
            
        probes.append( probe )

# =============================================================================
# plot
# =============================================================================

os.chdir( outpath )

if premultiply:
    if normalize: ylabel = r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$'
    else:         ylabel = r'$f \cdot PSD(f)$'
else:             ylabel = r'$PSD(f)$'

fig, ax = plt.subplots( figsize=(9, 8), constrained_layout=True )

for casenr in cases:
    
    probe = probes[casenr]
    
    if independent_len: lsep = len_sep[casenr]
    else:               lsep = len_sep[0]
    
    st = probe.psd_df['freq'] * 5.2 / 507 * lsep
    
    ax.semilogx( st, probe.psd_df[var],
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
ax.set_ylabel( ylabel )

ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(2)

if showlegend: ax.legend()

if premultiply: figname = 'pm_'    + figname
if normalize:   figname = 'norm_' + figname

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
    
    df = pd.DataFrame( { 'st': st, 'psd': probe.psd_df[var] } )
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
# ax.set_ylim( [0.0, 0.5] )

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
ax.set_ylabel( ylabel )

ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(2)

if showlegend: ax.legend()

# plt.savefig( figname + '_bin' + fmat )
plt.show()
plt.close()