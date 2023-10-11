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

snappath = datapath.split('/results')[0]+'/snapshots'
snapfiles = get_filelist( snappath, 'snapshot.bin' ) 

outpath  = datapath + outfolder
parametersfile = datapath.split('/results')[0] + '/case_parameters'


os.chdir(datapath)

spectra_files = get_filelist( snappath, 'spectra.dat' )

df_spectras = []
for spectra_file in spectra_files:
    
    df = pd.read_csv(spectra_file,sep='\s+')
    
    df_spectras.append( df )
    
    print(df)

E_uu = np.array([np.array(df['E_uu']) for df in df_spectras]).mean(axis=0)

print(E_uu)

fig,ax = plt.subplots()
ax.loglog(df['k_z'][:64],E_uu[:64])

# x = np.linspace(1,20)
# y = 1E6*x**(-5/3)

# ax.loglog(x,y)


    
plt.show()