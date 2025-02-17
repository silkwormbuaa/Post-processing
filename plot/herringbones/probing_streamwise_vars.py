#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   probing_streamwise_vars.py
@Time    :   2025/02/14 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import pickle
import numpy              as     np
import pandas             as     pd
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import LineData
plt.rcParams['font.size']           = 20

file0 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/wall_projection/streamwise_vars.pkl'
file  = '/home/wencan/temp/250120/postprocess/statistics/probing_x/probing_1.dat'

# =============================================================================

df = pd.read_csv( file, delimiter=r'\s+')

u_ref     = 507.0
p_ref     = 45447.289
rho_ref   = 0.9886
T_ref     = 160.15
delta_0   = 5.2

x = (np.array( df['x'] ) - 50.4)/delta_0

# =============================================================================

line = LineData()
with open( file0, 'rb' ) as f:    line.df = pickle.load( f )

# =============================================================================

p_fluc = np.sqrt(np.array(df['pp']) - np.array(df['p'])**2)/p_ref

fig, ax = plt.subplots( figsize=(15, 10) )

ax.plot( x, p_fluc, label='p_fluc', color='k', lw=2 )
ax.plot( line.df['x'], line.df['p_fluc'], label='smoothwall', color='r', lw=2 )

ax.set_xlabel(r'$(x-x_{\delta_0})/\delta_0$')
ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_\infty$")

ax.set_xlim(-20, 10)
# ax.set_ylim(0.0, 0.1)

plt.show()
plt.close()


