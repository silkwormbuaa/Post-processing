#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   analyse_pressure_signal.py
@Time    :   2025/02/04 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import numpy               as     np
import matplotlib.pyplot   as     plt
from   scipy.stats         import norm
from   matplotlib.gridspec import GridSpec

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.probe         import ProbeData
from   vista.timer         import timer
from   vista.params        import Params
from   vista.directories   import Directories
from   vista.tools         import get_filelist
from   vista.directories   import create_folder

# =============================================================================

casefolder = '/home/wencan/temp/smooth_mid'
vars       = ['p']
colors     = ['black']*len(vars)
prb_type   = 'pfmax'
time_range = [20.0,61.0]

# =============================================================================

dirs       = Directories( casefolder )
params     = Params( dirs.case_para_file )
prbfile    = dirs.fetch_prb_from_type( prb_type )

figname    = prb_type + '.png'

os.chdir( create_folder(dirs.pp_signal) )

# --- read in the probe data at the pressure fluctuation maximum point

prb_data   = ProbeData( prbfile, withT=params.prb_withT )
prb_data.cleandata( t_start=time_range[0], t_end=time_range[1] )

# --- compute statistics of signals

for var in vars:
    
    time    = np.array( prb_data.df['time'])
    signals = np.array( prb_data.df[var] )
    mean    = np.mean(  signals )
    std     = np.std(   signals )
    print(std)
    signals = (signals - mean)/std
    
    hist, bin_edges = np.histogram(signals, bins=100, density=True)
    
    bin_centers = (bin_edges[:-1]+bin_edges[1:])/2
    stdpdf      = norm.pdf( bin_centers, 0, 1)
    
    fig = plt.figure( figsize=(18,3) )
    
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)
    
    gs  = GridSpec( 1, 18 )
    
    ax1  = fig.add_subplot( gs[0,0:14])
    ax1.plot( time, signals, lw=0.5, color='black' )
    ax1.set_xlim( time_range )
    ax1.set_ylabel( var )
    
    ax2  = fig.add_subplot( gs[0,15:]) 
    ax2.plot( bin_centers, stdpdf, linestyle=':', color='gray')
    ax2.plot( bin_centers, hist,   linestyle='-', color='black')
    
    ax2.set_xlim([-3,3])
    ax2.set_ylim([0.0, 0.65])
    
    plt.savefig( var + figname )
    plt.show()
    plt.close()


