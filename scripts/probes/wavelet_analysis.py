
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   wavelet_analysis.py
@Time    :   2024/09/16 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   do wavelet analysis for the pressure signal
'''

import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt
from   scipy             import signal

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.directories import create_folder

# =============================================================================
# probe files

file = '/home/wencanwu/temp/231124_cwt_analysis/probe_00047.dat'
outpath = '/home/wencanwu/temp/231124_cwt_analysis/'

# =============================================================================

os.chdir( create_folder(outpath) )

prb = ProbeData( file, withT=True,  step=10 )
prb.cleandata( 20.0 )

time     = np.array( prb.df['time'] )
pressure = np.array( prb.df['p'] )

num_signal = len( pressure )
fs = num_signal/( time[-1] - time[0] )

print( num_signal )
print( fs )

wavelet = signal.morlet

width = np.arange(1,1000)

cwt_matrix = signal.cwt( pressure, wavelet, width )

# =============================================================================
# plot
# =============================================================================

plt.figure( figsize=(15, 10) )

plt.subplot(2,1,1)
plt.plot( time, pressure, 'k', lw=0.1 )
plt.xlim( time[0], time[-1] )
plt.title('Pressure signal')

plt.subplot(2,1,2)
plt.imshow( np.abs(cwt_matrix), aspect='auto', extent=[time[0], time[-1], width[0], width[-1]], cmap='coolwarm' )
plt.title('cwt of pressure signal')
plt.ylabel('width')
plt.xlabel('time/s')

plt.tight_layout()
plt.show()
plt.close()


