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

source_dir = os.path.dirname( os.path.dirname( os.path.realpath(__file__) ))
sys.path.append( source_dir )

import numpy             as np

import pandas            as pd

import matplotlib.pyplot as plt

from   pydmd             import DMD, SpDMD

from   pydmd.plotter     import plot_eigs

from   vista.timer       import timer

from   vista.tools       import get_filelist

from   vista.plot_style  import plot_eigens

from   vista.plot_style  import plot_amp_st

from   vista.snapshot    import Snapshot


snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots/snapshot_test/snapshots'

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

p_snap = np.array(p_snap).T/45447.289

pdmd = DMD(svd_rank=-1)

spdmd = SpDMD(svd_rank=-1,gamma=150000.0,rho=1.0)

with timer("DMD"):
    
    pdmd.fit(p_snap)
    
    spdmd.fit(p_snap)
    
pdmd.original_time['dt']=0.005
# plot eigenvalues

plot_eigens(pdmd.eigs)

Lsep = 58.9713176

velocity = 507.0

St = pdmd.frequency * Lsep / velocity

plot_amp_st( St, 
             np.abs(pdmd.amplitudes), 
             amp2 = np.abs(spdmd.amplitudes), 
             filename='serial.png' )

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
    
        
    
    
