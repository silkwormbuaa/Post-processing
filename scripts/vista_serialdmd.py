# -*- coding: utf-8 -*-
'''
@File    :   vista_dmd_serial.py
@Time    :   2023/05/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   computing DMD in serial
'''

import os

import sys

import numpy             as np

import pandas            as pd

import matplotlib.pyplot as plt

from   pydmd             import DMD

from   pydmd.plotter     import plot_eigs

sys.path.append('..')

from   utils.timer       import timer

from   utils.tools       import get_filelist

from   utils.plot_style  import plot_eigens

from   utils.plot_style  import plot_amp_st

from   vista.snapshot    import Snapshot


snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots/snapshot_01'

files = get_filelist( snap_dir )

# for file in files: print(file)

with timer("read in data "):
    
    p_snap = []
    
    x_snap = []
    
    z_snap = []

    for file in files:
        
        snapshot = Snapshot( file )
        
        snapshot.verbose = False
        
        snapshot.read_snapshot()
        
        snapshot.drop_ghost( buff=3 )
        
        
        # Assemble blocks data to one whole snapshot
        # (sorting is included in self.assemble_block() )
        
        snapshot.assemble_block()
        
        p = list( snapshot.df['p'] )
        
        p_snap.append( p )
        
        
    x_snap = np.array( snapshot.df['x'] )
    
    z_snap = np.array( snapshot.df['z'] )
    
    Nx     = np.size( np.unique(x_snap) )
    
    Nz     = np.size( np.unique(z_snap) )
    
    x_snap = x_snap.reshape( (Nz,Nx) ).T
    
    z_snap = x_snap.reshape( (Nz,Nx) ).T
        
print( "%d snapshots are loaded."%len(p_snap) )

p_snap = np.array(p_snap).T

pdmd = DMD(svd_rank=-1)



with timer("DMD"):
    
    pdmd.fit(p_snap)
    
pdmd.original_time['dt']=0.01
# plot eigenvalues

plot_eigens(pdmd.eigs)

Lsep = 50.069

velocity = 507.0

St = pdmd.frequency * Lsep / velocity

plot_amp_st( St, np.abs(pdmd.amplitudes) )

"""
# plot modes

for mode in pdmd.modes.T : 
    
    mode_real = np.array(mode.real).reshape( Nz, Nx ).T
    
    plt.pcolor(x_snap,z_snap,mode_real)
    
    plt.show()

t = np.linspace(0,1.44,145)


#for dynamic in pdmd.dynamics:
#    
#    plt.plot(t, dynamic.real)
#    plt.title("Dynamics")

print(type(pdmd.reconstructed_data))

print(np.shape(pdmd.reconstructed_data))

recon_data = pdmd.reconstructed_data.T

plt.pcolor(x_snap,z_snap,recon_data[-1].reshape(Nz,Nx).T.real)
    
plt.show()
    
"""
    
        
    
    
