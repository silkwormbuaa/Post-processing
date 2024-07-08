#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_psd_lines_free.py
@Time    :   2024/06/24 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Scripts for comparing arbitrary combination of psd lines
'''


import os
import sys
import pandas            as     pd
import numpy             as     np

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

import matplotlib.pyplot as     plt
from   vista.probe       import ProbeData
from   vista.timer       import timer
from   vista.line        import data_bin_avg

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 30

# =============================================================================

loc = 'sep'
independent_len = True
outpath         = '/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe/psd'
figname         = 'psd_lines_034'
fmt             = '.png'
cases           = [0,3,4]
showlegend      = False
normalize       = True

# =============================================================================

casescodes = ['smooth_adiabatic', '221014', '220926',
              '220825',           '220927', '221221',
              'smooth_isothermal']
sep_files  = ['probe_00144.dat', 'probe_00144.dat', 'probe_00145.dat',
              'probe_00135.dat',  'probe_00118.dat', 'probe_00175.dat',
              'probe_00144.dat']
pfmax_files= ['probe_00158.dat', 'probe_00152.dat', 'probe_00154.dat',
              'probe_00140.dat', 'probe_00125.dat', 'probe_00195.dat',
              'probe_00158.dat']
len_sep    = [9.52,               9.805522,       9.676337,
              11.340638,          13.12627,       13.266,
              9.52]
labels     = ['smooth',              r"$D/{\delta}_0=2.0$",  r"$D/{\delta}_0=1.0$",
               r"$D/{\delta}_0=0.5$",r"$D/{\delta}_0=0.25$", r"$D/{\delta}_0=0.125$",
               r"$\mathrm{isothermal}$"]
lstyles    = ['--',  ':',     '-.',    
              (0, (3, 1, 1, 1, 1, 1)), (0, (10, 3)), '-',     
              'dotted']
widths     = [4.0, 4.0, 4.0,
              4.0, 4.0, 4.0,
              4.0]
colors     = ['gray', 'green',  'blue', 
              'black', 'red',  'purple', 
              'yellow']
withTs     = [True,  False, False, 
              False, False, True, 
              False]
datapaths  = ['/media/wencanwu/Seagate Expansion Drive1/temp/'+casecode+'/probes/' 
             for casecode in casescodes]

# =============================================================================

if loc == 'sep':      prbfiles = sep_files
elif loc == 'pf_max': prbfiles = pfmax_files

probes = []

with timer("reading probes"):
    
    for casenr in cases:
        
        os.chdir( datapaths[casenr] )
        probe = ProbeData( prbfiles[casenr], withT=withTs[casenr] )
        probe.cleandata( t_start=20.0 )
        probe.get_fluc( 'p' )
        
        if normalize: 
            probe.pre_multi_psd( 'p_fluc', n_seg=8, overlap=0.5 )
            var = 'pmpsd_p_fluc'
        else:         
            probe.compute_psd( 'p_fluc', n_seg=8, overlap=0.5 )
            var = 'psd_p_fluc'
            
        probes.append( probe )

# =============================================================================
# plot
# =============================================================================

os.chdir( outpath )

fig, ax = plt.subplots( figsize=(12, 6), constrained_layout=True )

for casenr in cases:
    
    probe = probes[cases.index(casenr)]
    
    if independent_len: lsep = len_sep[casenr]
    else:               lsep = len_sep[0]
    
    st = probe.psd_df['freq'] * 5.2 / 507 * lsep
    
    ax.semilogx( st, probe.psd_df[var],
                 color=colors[casenr], 
                 linestyle=lstyles[casenr], 
                 linewidth=widths[casenr], 
                 label=labels[casenr] )
    
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

plt.savefig( figname + fmt )
plt.close()
    
# =============================================================================
# plot bin chart
# =============================================================================

fig, ax = plt.subplots( figsize=(12, 6), constrained_layout=True )

for casenr in cases:
    
    probe = probes[cases.index(casenr)]
    
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
             color=colors[casenr], 
             linewidth=widths[casenr], 
             linestyle=lstyles[casenr],
             label=labels[casenr] )

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
ax.set_ylabel( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$' )

ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(2)

if showlegend: ax.legend()

plt.savefig( figname + '_bin' + fmt )
plt.close()    
    

