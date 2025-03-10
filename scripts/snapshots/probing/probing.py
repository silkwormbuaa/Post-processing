#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   probing.py
@Time    :   2025/03/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   extract a probe line
'''


import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.directories import create_folder

# =============================================================================

case_dir   = '/home/wencan/temp/241030'
probe_type = 'X'
loc        = [ -1.99, 2.0 ]    # (z,y) in normalized coordinates
vars_in    = ['u', 'v', 'w', 'T', 'p']

# =============================================================================

dirs      = Directories( case_dir )

os.chdir( create_folder(dirs.pp_snp_probe) )

outfile   = f"probe_{probe_type}_{loc[0]:4.2f}_{loc[1]:4.2f}.dat"

# - read in case parameters

params    = Params( dirs.case_para_file )
delta     = params.delta_0

# - read in grid info

grid3d = GridData( dirs.grid )
grid3d.read_grid()

if   probe_type == 'X':
    loc = np.array(loc) * delta
elif probe_type == 'Y':
    loc[0] = loc[0] * delta + params.x_imp
    loc[1] = loc[1] * delta
elif probe_type == 'Z':
    loc[0] = loc[0] * delta + params.x_imp
    loc[1] = loc[1] * delta
    
bbox = params.snap_range_x + params.snap_range_y + params.snap_range_z 

print( f"bbox = {bbox}" )

blocklist, indx_probe = grid3d.select_probed_blockgrids( probe_type, loc, bbox )

print( f"blocklist = {blocklist}" )
print( f"indx_probe = {indx_probe}" )

# - find the snapshot

snappaths = dirs.snaps 
snap      = Snapshot( snappaths[0] + '/snapshot.bin' )
snap.read_snapshot( blocklist, vars_in )

df_snap = snap.get_probed_df( blocklist, grid3d, indx_probe, probe_type )
print( df_snap )

df_snap.to_string( outfile,
                   header = True,
                   index  = False,
                   float_format='%15.7f',
                   justify='left' )
print( f"Probe data written to {dirs.pp_snp_probe + '/' + outfile}." )

fig, ax = plt.subplots( figsize=(15, 10) )
ax.plot( (df_snap['x']-50.4)/5.2, df_snap['p']/params.p_ref, label='p', color='k', lw=2 )
ax.set_xlabel(r'$(x-x_{\delta_0})/\delta_0$')
ax.set_ylabel(r"$p/p_\infty$")

plt.show()
plt.close()

fig, ax = plt.subplots( figsize=(15, 10) )
ax.plot( (df_snap['x']-50.4)/5.2, df_snap['u']/params.u_ref, label='u', color='k', lw=2 )
ax.set_xlabel(r'$(x-x_{\delta_0})/\delta_0$')
ax.set_ylabel(r"$u/u_\infty$")
plt.show()
plt.close() 