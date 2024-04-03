#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run4.py
@Time    :   2024/03/12
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Scripts for comparing psd line plot
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

# Option zone
# =============================================================================

output_nr = 3              # 1,2,3
loc       = 'sep'          # 'sep' or 'pf_max'
pure      = False

# =============================================================================


datapath0 = '/media/wencanwu/Seagate Expansion Drive1/temp/smooth_adiabatic/probes'
datapath1 = '/media/wencanwu/Seagate Expansion Drive1/temp/240211/probes/'
datapath2 = '/media/wencanwu/Seagate Expansion Drive1/temp/220927/probes/'
datapath3 = '/media/wencanwu/Seagate Expansion Drive1/temp/240210/probes/'

outpath  = '/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe_ridge_height/psd'


# smooth wall
with timer("reading smooth wall probe"):
    os.chdir( datapath0 )
    if loc == 'sep':      probefile_s  = 'probe_00144.dat'
    elif loc == 'pf_max': probefile_s  = 'probe_00158.dat'
    Lsep_s       = 9.52     #13.12627403
    probe_s = ProbeData( probefile_s, withT=True )

if output_nr == 1:
    #240211
    with timer("reading one probe"):
        os.chdir( datapath1 )
        if loc == 'sep':      probefile_r = 'probe_00049.dat'
        elif loc == 'pf_max': probefile_r = 'probe_00058.dat'
        Lsep_r      = 11.44          #13.12627403
        probe_r  = ProbeData( probefile_r, withT=True )
        fig_name = 'varingLsep_240211_' + loc
        label_r  = r"$H/{\delta}_0=0.05$" 

elif output_nr == 2:
    #0927
    with timer("reading one probe"):
        os.chdir( datapath2 )
        if loc == 'sep':      probefile_r = 'probe_00118.dat'
        elif loc == 'pf_max': probefile_r = 'probe_00125.dat'
        Lsep_r      = 13.12627
        probe_r  = ProbeData( probefile_r, withT=False )
        fig_name = 'varingLsep_220927_' + loc
        label_r  = r"$H/{\delta}_0=0.10$"

elif output_nr == 3:
    #240210
    with timer("reading one probe"):
        os.chdir( datapath3 )
        if loc == 'sep':      probefile_r = 'probe_00020.dat'
        elif loc == 'pf_max': probefile_r = 'probe_00048.dat'
        Lsep_r      = 17.1           #13.12627403
        probe_r  = ProbeData( probefile_r, withT=True ) 
        fig_name = 'varingLsep_240210_' + loc
        label_r  = r"$H/{\delta}_0=0.20$"
    
# =============================================================================
# start plotting
# =============================================================================

os.chdir(outpath)
    
with timer('psd'):
    probe_s.cleandata( t_start=20 )
    probe_s.get_fluc( 'p' )
    probe_s.pre_multi_psd( 'p_fluc', n_seg=8, overlap=0.5 )
    
    probe_r.cleandata( t_start=20 )
    probe_r.get_fluc( 'p' )
    probe_r.pre_multi_psd( 'p_fluc', n_seg=8, overlap=0.5 )
    
    fig, ax = plt.subplots( figsize=[9,4], constrained_layout=True )
    
    St_Lsep_s = probe_s.psd_df['freq'] * 5.2 / 507 * Lsep_s
    St_Lsep_r = probe_r.psd_df['freq'] * 5.2 / 507 * Lsep_r
    
    ax.semilogx(St_Lsep_s, 
                probe_s.psd_df['pmpsd_p_fluc'],
                'red', 
                linewidth=1,
                label='smooth')
    
    ax.semilogx(St_Lsep_r, 
                probe_r.psd_df['pmpsd_p_fluc'],
                'blue', 
                linewidth=1,
                label=label_r)
    
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
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))

    ax.set_xlim( [0.01,1] )
    ax.set_ylim( [0.0,0.3] )
    
    ax.set_xlabel( r'$f L_{sep}/u_{\infty}$' )
    ax.set_ylabel( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$' )

    if not pure:
        ax.legend() 
    
#    ax.grid(visible=True, which='both',axis='both',color='gray',
#            linestyle='--',linewidth=0.5)
        
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        fig_name = fig_name + '_pure'

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)

    plt.savefig( fig_name+'zoomed_padded' )
#    plt.show()
    plt.close()

# =============================================================================
# plotting bin chart
# =============================================================================

with timer('bin chart'):
    
    dfs = pd.DataFrame( { 'St_Lsep_s':St_Lsep_s,
                          'psd': probe_s.psd_df['pmpsd_p_fluc'] } )
    
    dfr = pd.DataFrame( { 'St_Lsep_r':St_Lsep_r,
                          'psd': probe_r.psd_df['pmpsd_p_fluc'] } )
    
    bin_edges = np.logspace( -2, 2, 41, endpoint=True )
    
    dfsbin = data_bin_avg( dfs, 'St_Lsep_s', bin_edges )
    dfrbin = data_bin_avg( dfr, 'St_Lsep_r', bin_edges )
    
    # clear data
    
    dfsbin.dropna( inplace=True )
    dfrbin.dropna( inplace=True )
    
    # plot
    
    fig, ax = plt.subplots( figsize=[10,6], constrained_layout=True )
    
    barwidth = np.diff( bin_edges )
    
#    ax.bar( dfsbin['St_Lsep_s_mid'], dfsbin['psd'], width=barwidth, alpha=0.5, color='gray', label='smooth' )
#    ax.bar( dfrbin['St_Lsep_r_mid'], dfrbin['psd'], width=barwidth, alpha=0.5, color='blue', label=label_r )

#    spl = make_interp_spline( dfsbin['St_Lsep_s_mid'], dfsbin['psd'], k=3 )
#    ax.plot( np.logspace( -2,2,100 ), spl(np.logspace(-2,2,100 )), color='red', linewidth=2, label='smooth' )

#    spl = make_interp_spline( dfrbin['St_Lsep_r_mid'], dfrbin['psd'], k=3 )
#    ax.plot( np.logspace( -2,2,100 ), spl(np.logspace( -2,2,100 )), color='blue', linewidth=2, label=label_r )
    
    ax.plot( dfsbin['St_Lsep_s_mid'], dfsbin['psd'], color='red', linewidth=2, label='smooth' )
    ax.plot( dfrbin['St_Lsep_r_mid'], dfrbin['psd'], color='blue', linewidth=2, label=label_r )
    
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
    
    if not pure: ax.legend()
    
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)
    
    plt.savefig( fig_name+'_bin_line_nospline' )
    plt.show()
    plt.close()    
    
    
    
    