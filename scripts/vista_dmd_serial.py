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

from   vista.snapshot    import Snapshot


snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots/snap_test'

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

pdmd = DMD(svd_rank=145)


with timer("DMD"):
    
    pdmd.fit(p_snap)
    
# plot eigenvalues

plot_eigs(pdmd, show_axes=True, show_unit_circle=True)

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
    
    
        
    
    
