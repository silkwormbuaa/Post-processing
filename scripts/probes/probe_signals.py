#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   probe_signals.py
@Time    :   2024/09/24 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.timer       import timer
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.tools       import read_case_parameter
from   vista.directories import create_folder

# =============================================================================

casefolder = '/media/wencan/Expansion/temp/231124'
vars       = ['u','v','w','rho', 'p', 'T']
colors     = ['black']*len(vars)

# =============================================================================

dirs = Directories( casefolder )

params = read_case_parameter( dirs.case_para_file )
withT  = True if params.get('prb_withT').lower() == 'true' else False

# -- get all probe files
prb_files = get_filelist( dirs.prb_dir, 'probe_' )

os.chdir( create_folder( dirs.pp_signals ) )

# ----- first get the index of probes signals at snapshot time points.

prb_data = ProbeData( prb_files[0], withT=withT )
prb_data.cleandata( t_start=20.0 )
index    = prb_data.time_index( np.linspace(20.0, 61.0, 4101) )

clock    = timer("Plot signals from probes")
for file in prb_files:
    
    prb    = ProbeData( file, withT=withT )
    prb.cleandata( 20.0 )
    prb.df = prb.df.iloc[index]
    
    fig, axs = plt.subplots( len(vars),1, figsize=(10, 3*len(vars)), constrained_layout=True )
    
    for i, var in enumerate(vars):
        axs[i].plot( prb.df['time'], prb.df[var], colors[i], lw=0.1, label=var )
        axs[i].set_xlabel('time/s')
        axs[i].set_ylabel(var)
    
    figname = f'{os.path.splitext(os.path.basename(file))[0]}.png'
    plt.savefig( figname )
    print(f"Saved {figname}.")
    plt.close()

    progress = (i+1)/len(prb_files)
    print(f"{i+1}/{len(prb_files) } is done."+ clock.remainder(progress))
    print("---------------------\n")
    sys.stdout.flush()