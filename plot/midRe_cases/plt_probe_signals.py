#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_pressure_fluc.py
@Time    :   2024/09/16 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   plot pressure signal from probe data
'''

import os
import sys
import pandas            as     pd
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.directories import create_folder

# =============================================================================
# probe files

file1 = '/home/wencan/temp/231124/probes/probe_00047.dat'
file2 = '/home/wencan/temp/231124/probes/probe_00067.dat'
file3 = '/home/wencan/temp/231124/probes/probe_00221.dat'
file4 = '/home/wencan/temp/data_luis/midRe_streamwise/probe_at_separation/probe_00182.dat'

outpath = '/home/wencan/temp/231124/postprocess/probes/u_fluc'
# =============================================================================

os.chdir( create_folder(outpath) )

prb1 = ProbeData( file1, withT=True,  step=20 )
prb2 = ProbeData( file2, withT=True,  step=20 )
prb3 = ProbeData( file3, withT=True,  step=20 )
prb4 = ProbeData( file4, withT=False, step=20 )
prb1.cleandata( 20.0 )
prb2.cleandata( 20.0 )
prb3.cleandata( 20.0 )
prb4.cleandata( 20.0 )

fig, axs = plt.subplots( 3,1,figsize=(10, 15), constrained_layout=True )

axs[0].plot( prb1.df['time'], prb1.df['u']/507, 'red',   lw=0.1, label='flucmax' )
axs[1].plot( prb2.df['time'], prb2.df['u']/507, 'blue',  lw=0.1, label='sep' )
axs[2].plot( prb3.df['time'], prb3.df['u']/507, 'black', lw=0.1, label='reatt' )
#ax.plot( prb4.df['time'], prb4.df['u']/507, 'green', lw=0.1, label='sep_smooth' )

for ax in axs:
    ax.set_xlabel('time/s')
    ax.set_ylabel(r'$u/u_{\infty}$')

plt.savefig( 'velocity_fluc_smooth.png' )
plt.show()
plt.close()