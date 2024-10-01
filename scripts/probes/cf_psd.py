#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   cf_psd.py
@Time    :   2024/09/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Analyse a signal probe's signal
'''

import os
import sys
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.directories import Directories
from   vista.line        import data_bin_avg
from   vista.tools       import read_case_parameter
from   vista.material    import get_visc
from   vista.directories import create_folder

# =============================================================================

dirpath = '/media/wencan/Expansion/temp/231124'
file    = dirpath + '/probes/probe_00070.dat'
outpath = dirpath + '/postprocess/cf_psd/'

# =============================================================================

dirs = Directories( dirpath )
os.chdir( create_folder(outpath) )

params   = read_case_parameter( dirs.case_para_file )
u_ref    = float(params.get('u_ref'))
rho_ref  = float(params.get('rho_ref'))
delta    = float(params.get('delta_0'))
lsep     = float(params.get('Lsep'))
Re_ref   = float(params.get('Re_ref'))
visc_law = params.get('visc_law')
p_dyn    = 0.5 * rho_ref * u_ref**2

prb = ProbeData( file, withT=True )
prb.cleandata( 20.0 )

walldist = prb.xyz[1]

time    = np.array( prb.df['time'] )
ts      = np.array( prb.df['T'] )
u       = np.array( prb.df['u'] )
mu      = get_visc( ts, Re_ref, law= visc_law )
cf      = mu*u/walldist/p_dyn

prb.add_data( 'cf', cf )
prb.get_fluc( ['cf'] )
prb.pre_multi_psd( 'cf_fluc',n_seg=8, overlap=0.5 )

# fig, ax = plt.subplots( figsize=(15, 10) )

# ax.plot( time[::10], cf[::10], 'k', lw=0.1 )
# print( np.mean(cf) )
# plt.show()
# plt.close()

fig, ax = plt.subplots( figsize=(15, 10) )

st = prb.psd_df['freq'] / u_ref * lsep

df = pd.DataFrame( {'st':st, 'cf':prb.psd_df['pmpsd_cf_fluc']} )
bin_edges = np.logspace( -2, 2, 41, endpoint=True ) 

df_bin = data_bin_avg(df, 'st', bin_edges )

df_bin.dropna( inplace=True )

ax.plot( df_bin['st'], df_bin['cf'], 'k', lw=2 )
ax.set_xscale( 'log' )
ax.set_xlim( [0.01,100] )

plt.show()