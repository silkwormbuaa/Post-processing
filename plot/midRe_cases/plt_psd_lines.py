#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_psd_4lines.py
@Time    :   2024/10/08 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import pandas            as     pd
import numpy             as     np
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0] 
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.directories import Directories
from   vista.line        import data_bin_avg
from   vista.directories import create_folder

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 30

# =============================================================================
# option zone
# =============================================================================

independent_len = False
figname         = 'psdlines_pfmax'
fmat            = '.png'
figsize         = (10,7)
cases_nr        = [0,1,2,3,4]
showlegend      = False
premultiply     = True
normalize       = False
outpath         = '/home/wencan/temp/DataPost/midRe/psd/jfm_second_revision/'

cases    = ['smooth_adiabatic','220927','smooth_mid','231124','241030']
color    = ['gray',             'orangered',       'black',           'steelblue',        'yellowgreen'          ,'red'] 
label    = [r'$\mathcal{LS}$',  r'$\mathcal{LR}$', r'$\mathcal{HS}$', r'$\mathcal{HR}1$', r'$\mathcal{HR}2$'     ,r'$\mathcal{HR}3$']
lstyle   = ['-',                '-.',               '--',             ':',                (0, (3, 1, 1, 1, 1, 1)),'--']
lwidth   = [4.0,              4.0,               4.0,             4.0,                4.0,                  4.0]
locs     = ['pfmax','pfmax','pfmax','pfmax','pfmax']

os.chdir( create_folder(outpath) )

# =============================================================================
# un-binned plot
# =============================================================================

fig, ax = plt.subplots( figsize=figsize, constrained_layout=True )

for nr in cases_nr:
    
    case     = cases[nr]
    casepath = '/home/wencan/temp/' + case
    dirs     = Directories( casepath )
    params   = Params( dirs.case_para_file )
    
    prbfile = dirs.fetch_prb_from_type( locs[nr] )

    with timer("reading probes"):
    
        probe = ProbeData( prbfile, withT=params.prb_withT )
        probe.cleandata( t_start=20.0 )
        probe.get_fluc( 'p' )
        
        if premultiply: 
            probe.pre_multi_psd( 'p_fluc', n_seg=8, overlap=0.5, normalize=normalize, nfft=None )
            var = 'pmpsd_p_fluc'
        else:         
            probe.compute_psd( 'p_fluc', n_seg=8, overlap=0.5 )
            var = 'psd_p_fluc'
            

    if premultiply:
        if normalize: ylabel = r'$f \cdot \mathcal{P}(f)/ \int \mathcal{P}(f) \mathrm{d} f$'
        else:         ylabel = r'$f \cdot \mathcal{P}(f)$'
    else:             ylabel = r'$\mathcal{P}(f)$'


    if independent_len: lsep = params.lsep
    else:               lsep = 9.52
    
    st = probe.psd_df['freq'] * 5.2 / 507 * lsep
    
    ax.semilogx( st, probe.psd_df[var],
                 color    =color[nr], 
                 linestyle=lstyle[nr], 
                 linewidth=lwidth[nr], 
                 label=label[nr] )
    
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


fig, ax = plt.subplots( figsize=figsize, constrained_layout=True )

for nr in cases_nr:
    
    case     = cases[nr]
    casepath = '/home/wencan/temp/' + case
    dirs     = Directories( casepath )
    params   = Params( dirs.case_para_file )
    
    prbfile = dirs.fetch_prb_from_type( locs[nr] )

    with timer("reading probes"):
    
        probe = ProbeData( prbfile, withT=params.prb_withT )
        probe.cleandata( t_start=20.0 )
        probe.get_fluc( 'p' )
        
        if premultiply: 
            probe.pre_multi_psd( 'p_fluc', n_seg=8, overlap=0.5, normalize=normalize )
            var = 'pmpsd_p_fluc'
        else:         
            probe.compute_psd( 'p_fluc', n_seg=8, overlap=0.5 )
            var = 'psd_p_fluc'
            

    if premultiply:
        if normalize: ylabel = r'$f \cdot \mathcal{P}(f)/ \int \mathcal{P}(f) \mathrm{d} f$'
        else:         ylabel = r'$f \cdot \mathcal{P}(f)$'
    else:             ylabel = r'$\mathcal{P}(f)$'


    if independent_len: lsep = params.lsep
    else:               lsep = 9.52
    
    st = probe.psd_df['freq'] * 5.2 / 507 * lsep
    
    df = pd.DataFrame( { 'st': st, 'psd': probe.psd_df[var] } )
    bin_edges = np.logspace( -2, 2, 41, endpoint=True )
    
    df_bin = data_bin_avg( df, 'st', bin_edges )
    
    # clear data
    
    df_bin.dropna( inplace=True )
    
    # plot    
    
    barwidth = np.diff( bin_edges )
    
    ax.plot( df_bin['st_mid'], df_bin['psd'], 
             color=color[nr], 
             linewidth=lwidth[nr], 
             linestyle=lstyle[nr],
             label=label[nr] )

ax.vlines( 0.4, -1.0, 8.0e6, linestyles='--', linewidth=2, color='black' )

ax.set_xscale( 'log' )
ax.set_xlim( [0.01,100] )
if normalize: 
    ax.set_ylim( [0.0,0.5] )
else:
    ax.set_ylim( [0.0,8.0e6] )

ax.set_position([0.10,0.18,0.85,0.75])
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

ax.set_xlabel( r'$St_{Lsep} = f L_{sep}/u_{\infty}$' )
ax.set_ylabel( ylabel )

ax.spines[:].set_color('black')
ax.spines[:].set_linewidth(2)

if showlegend: ax.legend()

plt.savefig( figname + '_bin' + fmat )
plt.close()