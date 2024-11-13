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
from   scipy.stats       import skew, kurtosis

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.params      import Params
from   vista.directories import Directories
from   vista.directories import create_folder

# =============================================================================
# probe files

outpath = '/home/wencan/temp/231124/postprocess/probes/u_fluc'

file1 = '/home/wencan/temp/231124/probes/probe_00047.dat'
file2 = '/home/wencan/temp/231124/probes/probe_00077.dat'  #67
file3 = '/home/wencan/temp/231124/probes/probe_00221.dat'

casedir = '/home/wencan/temp/231124'
dirs = Directories( casedir )
pfmaxfile = dirs.fetch_prb_from_type( 'pfmax' )

prb_file_list = [ file1, file2, file3 ]

prb_list = list()

# =============================================================================

os.chdir( create_folder(outpath) )

for prb_file in prb_file_list:
    prb = ProbeData( prb_file, withT=True )
    prb.cleandata( 20.0 )
    prb_list.append( prb )

fig, axs = plt.subplots( 3,1,figsize=(10, 15), constrained_layout=True )

for i, ax in enumerate( axs ):
    ax.plot( prb_list[i].df['time']/507, prb_list[i].df['u'], 'black', lw=0.1)
    ax.set_xlabel('time/s')
    ax.set_ylabel(r'$u/u_{\infty}$')

plt.savefig( 'velocity_fluc_smooth.png' )
plt.show()
print( skew(prb_list[0].df['u']), kurtosis(prb_list[0].df['u']) )
print( skew(prb_list[1].df['u']), kurtosis(prb_list[1].df['u']) )
print( skew(prb_list[2].df['u']), kurtosis(prb_list[2].df['u']) )