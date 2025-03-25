#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   shock_cross_correlation.py
@Time    :   2024/11/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Cross-correlation analysis of the shock motion at different wall 
             normal locations.
             Now this script can be replaced by ../correlation_analysis/cross_correlation_analysis.py
'''

import os
import sys
import numpy              as     np
import pandas             as     pd
import matplotlib.pyplot  as     plt
from   scipy.signal       import correlate

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(fontsize=15)

# =============================================================================

casedir = '/home/wencan/temp/231124'

shockmotion1 = casedir + '/postprocess/snapshots/shock_tracking/2d_y2/shock_tracking2d.csv'
shockmotion2 = casedir + '/postprocess/snapshots/shock_tracking/2d_y6/shock_tracking2d.csv'

# =============================================================================

df1 = pd.read_csv( shockmotion1, delimiter=r'\s+' )
df2 = pd.read_csv( shockmotion2, delimiter=r'\s+' )

print(df1['x_shock'])
print(df2['x_shock'])

x1_fluc = df1['x_shock'] - df1['x_shock'].mean()
x2_fluc = df2['x_shock'] - df2['x_shock'].mean()

correlation = correlate( x1_fluc, x2_fluc, mode='full', method='auto' )
correlation /= (np.std(x1_fluc) * np.std(x2_fluc) * len(x1_fluc))

lags = np.arange(-len(df1['x_shock'])+1, len(df1['x_shock']))*0.01

lag_at_max_corr = lags[np.argmax(correlation)]

print(f"lag at max correlation: {lag_at_max_corr:.2f} ms")

# =============================================================================
# plot

fig, ax = plt.subplots( figsize=(12,6) )
ax.plot( lags, correlation )
ax.axvline( lag_at_max_corr, color='r', ls='--' )
plt.show()
plt.close()
