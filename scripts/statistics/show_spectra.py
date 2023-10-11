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

from   vista.statistic   import StatisticData

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.tools       import read_case_parameter
from   vista.tools       import get_filelist

from   vista.line        import LineData


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

parameters = read_case_parameter( parametersfile )
delta   = float( parameters.get('delta_0') )
u_ref   = float( parameters.get('u_ref') )


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


# -- Euu
fig,ax = plt.subplots()
ax.loglog( df['k_z'][1:64]*delta, E_uu[1:64]/(delta*u_ref**2) )

x = np.linspace(10,100,100)
y = 10*x**(-5/3)
ax.loglog(x,y)
    
plt.show()
plt.close()

# -- Evv
fig,ax = plt.subplots()
ax.loglog( df['k_z'][1:64]*delta, E_vv[1:64]/(delta*u_ref**2) )

x = np.linspace(20,100,100)
y = 3*x**(-5/3)
ax.loglog(x,y)

plt.show()
plt.close()

# -- Eww
fig,ax = plt.subplots()
ax.loglog( df['k_z'][1:64]*delta, E_ww[1:64]/(delta*u_ref**2) )

x = np.linspace(10,100,100)
y = 3*x**(-5/3)
ax.loglog(x,y)

plt.show()
plt.close()

# -- Ek
fig,ax = plt.subplots()
ax.loglog( df['k_z'][1:64]*delta, E_k[1:64]/(delta*u_ref**2) )

x = np.linspace(8,100,100)
y = 25*x**(-5/3)
ax.loglog(x,y)

plt.show()
plt.close()