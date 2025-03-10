#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_spectra.py
@Time    :   2023/10/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import pandas            as     pd
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.tools       import get_filelist
from   vista.params      import Params

# =============================================================================

outfolder = '/spectra'
probe_type = 'Z'
loc = (-20.0,0.52)

# =============================================================================

datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'

spectrapath = datapath.split('/results')[0]+'/spectra'
snapfiles = get_filelist( spectrapath, 'snapshot.bin' ) 

outpath  = datapath + outfolder
parametersfile = datapath.split('/results')[0] + '/case_parameters'

parameters = Params( parametersfile )
delta      = parameters.delta_0
u_ref      = parameters.u_ref


os.chdir(datapath)

spectra_files = get_filelist( spectrapath, 'spectra.dat' )

df_spectras = []
for spectra_file in spectra_files:
    
    df = pd.read_csv(spectra_file,sep='\s+')
    
    df_spectras.append( df )
    
    print(f"read in spectra of {spectra_file}")

E_uu = np.array([np.array(df['E_uu']) for df in df_spectras]).mean(axis=0)
E_vv = np.array([np.array(df['E_vv']) for df in df_spectras]).mean(axis=0)
E_ww = np.array([np.array(df['E_ww']) for df in df_spectras]).mean(axis=0)

E_k = E_uu + E_vv + E_ww

print(df['k_z'])

os.chdir( datapath+outfolder )

# -- Euu
fig,ax = plt.subplots(figsize=(5,10))
ax.loglog( df['k_z'][1:64]*delta, E_uu[1:64]/(delta*u_ref**2) )

x = np.linspace(10,100,100)
y = 0.5*x**(-5/3)
ax.loglog(x,y,ls='--')
    
plt.savefig("Euu")
plt.close()

# -- Evv
fig,ax = plt.subplots(figsize=(5,10))
ax.loglog( df['k_z'][1:64]*delta, E_vv[1:64]/(delta*u_ref**2) )

x = np.linspace(20,100,100)
y = 0.3*x**(-5/3)
ax.loglog(x,y,ls='--')

plt.savefig("Evv")
plt.close()

# -- Eww
fig,ax = plt.subplots(figsize=(5,10))
ax.loglog( df['k_z'][1:64]*delta, E_ww[1:64]/(delta*u_ref**2) )

x = np.linspace(10,100,100)
y = 0.3*x**(-5/3)
ax.loglog(x,y,ls='--')

plt.savefig("Eww")
plt.close()

# -- Ek
fig,ax = plt.subplots(figsize=(5,10))
ax.loglog( df['k_z'][1:64]*delta, E_k[1:64]/(delta*u_ref**2) )

x = np.linspace(8,100,100)
y = 0.6*x**(-5/3)
ax.loglog(x,y,ls='--')

plt.savefig("Ek")
plt.close()