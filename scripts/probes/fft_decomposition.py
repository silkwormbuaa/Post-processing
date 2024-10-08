#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   fft_decomposition.py
@Time    :   2024/10/07 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Using fft to decompose different frequency components of pressure signals
'''

import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData


# =============================================================================

probefile = '/home/wencanwu/test/231124/probes/probe_00047.dat'

#probefile = '/home/wencanwu/test/smooth_mid/probes/probe_00084.dat'

# =============================================================================

prb = ProbeData( probefile, withT=True, step=20 )
prb.cleandata( 20.0 )


time = np.array( prb.df['time'] )
p    = np.array( prb.df['p'] )
p    = p - np.mean(p)

print( len(time) )

# -- fft

tf   = np.fft.fftfreq(len(time), d=time[1]-time[0])
pf   = np.fft.fft(p)

# -- inverse fft

f_cutoff = 3.0
p1   = pf.copy()
p1[np.abs(tf) > f_cutoff] = 0.0
pl_re = np.fft.ifft(p1)

p2   = pf.copy()
p2[np.abs(tf) < f_cutoff] = 0.0
ph_re = np.fft.ifft(p2)

fig, ax = plt.subplots( 3,1, figsize=(10, 10))

ax[0].plot( time, p/45447.289, lw=0.2, label='p' )
ax[0].plot( time, pl_re/45447.289, 'r', lw=0.5)
ax[0].set_ylim([-0.3,0.3])
ax[0].set_xlabel('Time [ms]')
ax[0].set_ylabel(r'$p/p_{\infty}$')
ax[0].grid()

ax[1].plot( tf[:len(time)//2]*48/507, np.abs(pf[:len(time)//2])/(len(pf)//2)/45447.289, lw=0.2, label='p' )
ax[1].set_ylim([-0.001,0.04])
ax[1].set_xscale('log')
ax[1].set_xlabel(r'$St_{sep}$')
ax[1].set_ylabel('Amplitude')

ax[2].plot( time, ph_re/45447.289, 'b', lw=0.2)
ax[2].set_ylim([-0.25,0.25])

print( np.std(pl_re), np.std(ph_re) )

plt.show()
