#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   cross_corelation_analysis.py
@Time    :   2024/09/11
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   joint analysis of shock motion, bubble size and instantaneous pressure
'''

import os
import sys
import pickle
import numpy              as     np
import pandas             as     pd
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.timer        import timer
from   vista.params       import Params
from   vista.directories  import Directories
from   vista.directories  import create_folder
from   vista.tools        import get_filelist
from   vista.tools        import create_linear_interpolator

#set_plt_rcparams()
# =============================================================================

case_dir = '/media/wencan/Expansion/temp/241030'

shockpath3d  = case_dir + '/postprocess/shock/group1'
bubblepath   = case_dir + '/postprocess/bubble/bubble_size.dat'
pressurepath = case_dir + '/postprocess/probes/pre_ridge/pressure_ridge.pkl'

outputpath   = case_dir + '/postprocess/joint_analysis'
os.chdir( create_folder( outputpath ) )

dirs = Directories( case_dir )
# =============================================================================

# -----------------------------------------------------------------------------
# - read in case parameters
# -----------------------------------------------------------------------------
params  = Params( dirs.case_para_file )
p_ref   = params.p_ref
h       = params.H
x_imp   = params.x_imp
delta_0 = params.delta_0
# the following parameters are normalized by delta_0
x_sep   = params.x_sep
x_att   = params.x_att
x_pfmax = params.x_pfmax

# -----------------------------------------------------------------------------
# - read in the shock motion data
# -----------------------------------------------------------------------------

shock3d_file = get_filelist( shockpath3d, 'shock_tracking1' )[0]

with open(shock3d_file, 'rb') as f:
    times      = pickle.load(f) 
    shocklines = pickle.load(f)

x_shocks     = list()
x_shocks_mid = list()

for shockline in shocklines:
    x_shock = np.array( shockline['x'] )
    x_shocks.append( x_shock.mean() )
    x_shocks_mid.append( x_shock[int(len(x_shock)/2)] )
    
x_mean = np.mean(x_shocks)
x_fluc = np.array(x_shocks - x_mean)

# -----------------------------------------------------------------------------
# - read in the bubble size data ['itstep', 'itime', 'bubble_volume']
# -----------------------------------------------------------------------------

bubble_df = pd.read_csv(bubblepath,delimiter=r'\s+')

bb_size = np.array(bubble_df['bubble_volume'])
bb_size_mean = bb_size.mean()
bb_size_fluc = bb_size - bb_size_mean

times = ( np.array(times) - 20.0 ) * 507.0 * (1.0/7.15)

norm_shock  = 0.5*(np.array(x_fluc).max() - np.array(x_fluc).min())
norm_bbsize = 0.5*(np.array(bb_size_fluc).max() - np.array(bb_size_fluc).min())

# -----------------------------------------------------------------------------
# - read in the streamwise pressure data
# -----------------------------------------------------------------------------

with open(pressurepath, 'rb') as f:
    times, x_locs, pres = pickle.load( f )

pres = np.array( pres )
n_locs, n_time = pres.shape
x_locs = (np.array(x_locs) - 50.4) / 5.2

# - exclude two probes accidentially placed at the ridge
pres   = pres[:-2,:]
x_locs = x_locs[:-2]

# -----------------------------------------------------------------------------
# - read in the mean pressure data
# -----------------------------------------------------------------------------

df = pickle.load( open(dirs.pp_wall_proj+'/streamwise_vars.pkl', 'rb') )

# -----------------------------------------------------------------------------
# - plot bubble size, shock motion and instantenous pressure
# -----------------------------------------------------------------------------

# - plot shock motion

fin = create_linear_interpolator( df['x'], df['Cp'] )
clock = timer('plotting')

for i in range( len(times) ): 

    fig, axs = plt.subplots( 3, 1, figsize=(20,20) )
    
    axs[0].plot( times, x_fluc/norm_shock,'b' )
    axs[1].plot( times, bb_size_fluc/norm_bbsize,'r', ls=':' )

    #ax.plot( times, np.array(x_shocks_mid)/7.15,'r', ls=':' )
    axs[0].set_xlim( 20, 61 )
    axs[0].set_ylim( -1.2, 1.2 )
    axs[0].set_title('Spanwise averaged shock location')
    axs[1].set_xlim( 20, 61 )
    axs[1].set_ylim( -1.2, 1.2 )
    axs[1].set_title('bubble size')

    # add dot of the current time
    
    axs[0].plot( times[i], x_fluc[i]/norm_shock, 'ro' )
    axs[1].plot( times[i], bb_size_fluc[i]/norm_bbsize, 'ro' )

    # - plot pressure distribution
    axs[2].plot( df['x'], df['Cp'], color='red' )
    axs[2].plot( x_locs, pres[:,i]/p_ref, '--' )
    axs[2].plot([x_sep, x_att],[fin(x_sep), fin(x_att)], color='blue', marker='p',markersize=10,ls='')
    axs[2].plot([x_pfmax],[fin(x_pfmax)], color='red', marker='p',markersize=10,ls='')
    axs[2].set_title('Instantaneous pressure distribution')
    axs[2].set_xlim( -15, 10 )
    axs[2].set_ylim( 0.9, 2.5 )
    
    fig.suptitle( f"t={times[i]:.2f}s" )
    
    plt.savefig( f"joint_analysis_{i:06d}.png" ) 
    
    print(f"{i+1}/{len(times)} is done. " + clock.remainder((i+1)/len(times)))  
#    plt.show()
    plt.close()
