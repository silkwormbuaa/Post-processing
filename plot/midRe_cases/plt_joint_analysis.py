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

from   vista.tools        import get_filelist
from   vista.plot_setting import set_plt_rcparams

#set_plt_rcparams()
# =============================================================================

shockpath3d = '/media/wencan/Expansion/temp/231124/postprocess/shock_tracking/3d'
bubblepath  = '/media/wencan/Expansion/temp/231124/postprocess/bubble_size/bubble_size.dat'

# =============================================================================

# - read in the shock motion data
# -----------------------------------------------------------------------------

shock3d_files = get_filelist( shockpath3d, 'shock_tracking3d' )

times = list()
shocklines = list()

for shock3d_file in shock3d_files:
    with open(shock3d_file, 'rb') as f:
        time = pickle.load(f) 
        shockline = pickle.load(f)

    times = times + time
    shocklines = shocklines + shockline

x_shocks = list()
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

times_bb = np.array(bubble_df['itime'])

bb_size = np.array(bubble_df['bubble_volume'])
bb_size_mean = bb_size.mean()
bb_size_fluc = bb_size - bb_size_mean

# - plot shock motion

fig, ax = plt.subplots( figsize=(15,8) )

times = ( np.array(times) - 20.0 ) * 507.0 * (1.0/7.15)

norm_shock = 0.5*(np.array(x_fluc).max() - np.array(x_fluc).min())
ax.plot( times, np.array(x_fluc)/norm_shock,'b' )

norm_bbsize = 0.5*(np.array(bb_size_fluc).max() - np.array(bb_size_fluc).min())
ax.plot( times, bb_size_fluc/norm_bbsize,'r', ls=':' )

#ax.plot( times, np.array(x_shocks_mid)/7.15,'r', ls=':' )
ax.set_ylim( -1.1, 1.1 )
ax.set_title('Spanwise averaged shock location')

plt.show()
plt.close()
