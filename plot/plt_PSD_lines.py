#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run4.py
@Time    :   2022/11/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Scripts for comparing psd line plot
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

from   vista.psd         import ProbeData

from   vista.timer       import timer

from   vista.tools       import get_filelist
from   vista.tools       import read_case_parameter

from   vista.line        import data_bin_avg

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 30

# Option zone
# =============================================================================

output_nr = 3              # 1,2,3,4,5
loc       = 'sep'          # 'sep' or 'pf_max'
pure      = False

# =============================================================================


datapath0 = '/home/wencanwu/my_simulation/temp/smooth_wall/probes/probe_x/'
datapath1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/probe_x/'
datapath2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/probes/probe_x/'
datapath3 = '/media/wencanwu/Seagate Expansion Drive1/temp/220825/probes/probe_x/'
datapath4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/probes/probe_x/'
datapath5 = '/media/wencanwu/Seagate Expansion Drive1/temp/221221/probes/probe_x/'

outpath  = '/home/wencanwu/my_simulation/temp/DataPost/PSD'


# smooth wall
with timer("reading smooth wall probe"):
    os.chdir( datapath0 )
    if loc == 'sep':      probefile_s  = 'probe_00142.dat'
    elif loc == 'pf_max': probefile_s  = 'probe_00156.dat'
    Lsep_s       = 13.12627403
    probe_s = ProbeData( probefile_s, withT=False )

if output_nr == 1:
    #1014
    with timer("reading one probe"):
        os.chdir( datapath1 )
        if loc == 'sep':      probefile_r = 'probe_00144.dat'
        elif loc == 'pf_max': probefile_r = 'probe_00152.dat'
        Lsep_r      = 9.805522
        probe_r  = ProbeData( probefile_r, withT=False )
        fig_name = 'psd_line_1_1014_' + loc
        label_r  = r"$D/{\delta}_0=2.0$" 

elif output_nr == 2:
    #0926
    with timer("reading one probe"):
        os.chdir( datapath2 )
        if loc == 'sep':      probefile_r = 'probe_00145.dat'
        elif loc == 'pf_max': probefile_r = 'probe_00154.dat'
        Lsep_r      = 9.676337
        probe_r  = ProbeData( probefile_r, withT=False )
        fig_name = 'psd_line_2_0926_' + loc
        label_r  = r"$D/{\delta}_0=1.0$" 

elif output_nr == 3:
    #0825
    with timer("reading one probe"):
        os.chdir( datapath3 )
        if loc == 'sep':      probefile_r = 'probe_00135.dat'
        elif loc == 'pf_max': probefile_r = 'probe_00140.dat'
        Lsep_r      = 11.340638
        probe_r  = ProbeData( probefile_r, withT=False )
        fig_name = 'psd_line_3_0825_' + loc
        label_r  = r"$D/{\delta}_0=0.5$" 

elif output_nr == 4:
    #0927
    with timer("reading one probe"):
        os.chdir( datapath4 )
        if loc == 'sep':      probefile_r = 'probe_00118.dat'
        elif loc == 'pf_max': probefile_r = 'probe_00125.dat'
        Lsep_r      = 13.12627
        probe_r  = ProbeData( probefile_r, withT=False )
        fig_name = 'psd_line_4_0927_' + loc
        label_r  = r"$D/{\delta}_0=0.25$"

elif output_nr == 5:
    #1221
    with timer("reading one probe"):
        os.chdir( datapath5 )
        if loc == 'sep':      probefile_r = 'probe_00175.dat'
        elif loc == 'pf_max': probefile_r = 'probe_00195.dat'
        Lsep_r      = 13.266
        probe_r  = ProbeData( probefile_r ) 
        fig_name = 'psd_line_5_1221_' + loc
        label_r  = r"$D/{\delta}_0=0.125$"


""" #1125
with timer("reading one probe"):
    os.chdir( folderpath1 )
#    probefile_r = 'probe_00176.dat'
    probefile_r = 'probe_00189.dat'
    Lsep_r      = 13.12627
    probe_r = ProbeData( probefile_r ) """
    
# =============================================================================
# start plotting
# =============================================================================

with timer("clean data"):
    probe_s.cleandata( starttime=20 )
    probe_r.cleandata( starttime=20 )
    
with timer('psd'):
    
    os.chdir(outpath)
    
    probe_s.psd( 8, 0.5 )
    probe_r.psd( 8, 0.5 )
    
    fig, ax = plt.subplots( figsize=[9,4], constrained_layout=True )
    
    St_Lsep_s = probe_s.St * Lsep_s
    St_Lsep_r = probe_r.St * Lsep_r
    
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

    ax.semilogx(St_Lsep_s, 
                probe_s.nd_pprime_fwpsd,
                'red', 
                linewidth=1,
                label='smooth')
    
    ax.semilogx(St_Lsep_r, 
                probe_r.nd_pprime_fwpsd,
                'blue', 
                linewidth=1,
                label=label_r)
    
    ax.set_xlim( [0.01,1] )
    ax.set_ylim( [0.0,0.3] )
    
#    ax.set_xlabel( r'$St_{sep}$' )
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
                          'psd': probe_s.nd_pprime_fwpsd } )
    
    dfr = pd.DataFrame( { 'St_Lsep_r':St_Lsep_r,
                          'psd': probe_r.nd_pprime_fwpsd } )
    
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

    spl = make_interp_spline( dfsbin['St_Lsep_s_mid'], dfsbin['psd'], k=3 )
    ax.plot( np.logspace( -2,2,100 ), spl(np.logspace(-2,2,100 )), color='red', linewidth=2, label='smooth' )

    spl = make_interp_spline( dfrbin['St_Lsep_r_mid'], dfrbin['psd'], k=3 )
    ax.plot( np.logspace( -2,2,100 ), spl(np.logspace( -2,2,100 )), color='blue', linewidth=2, label=label_r )
    
    ax.set_xscale( 'log' )
    ax.set_xlim( [0.01,100] )
    ax.set_ylim( [0.0, 0.4] )

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
    
    plt.savefig( fig_name+'_bin_line' )
    plt.show()
    plt.close()    
    
    
    
    