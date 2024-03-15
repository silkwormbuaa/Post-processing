#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run4.py
@Time    :   2022/11/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Scripts for getting PSD
'''

import os
import sys
import time

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.directories import Directories
from   vista.probe       import ProbeFile
from   vista.timer       import timer
from   vista.colors      import colors as col
from   vista.tools       import read_case_parameter


dirs   = Directories( os.getcwd() )
probes = ProbeFile( dirs.set_prb )

# -- enter working directory and read relevant parameters 

os.chdir( dirs.prb_dir )

prb_data = dirs.probes
n_data = len( prb_data )
print(f"We are in directory:{dirs.prb_dir}\n")
print(f"We have got {n_data:5d} probes data.\n")

parameters = read_case_parameter( dirs.case_para_file )
h = float( parameters.get('H') )
prb_withT = True if parameters.get('prb_withT').lower() == 'true' else False

# -- check if number of probe data and number of probes are consistent

if n_data != len( probes.probes ):
    raise ValueError("Number of probe data and number of probes are not consistent.")

# -- check if probe location and classify them to either ridge/valley or others

if not os.path.exists( dirs.pp_psd_ridge  ): os.mkdir( dirs.pp_psd_ridge  )
if not os.path.exists( dirs.pp_psd_valley ): os.mkdir( dirs.pp_psd_valley )
if not os.path.exists( dirs.pp_psd_others ): os.mkdir( dirs.pp_psd_others )

# -- start computing PSD 

with timer('Computing PSD '):
    
    # start 
    
    n_ridge = 0 ; n_valley = 0 ; n_others = 0
    
    for i, probe in enumerate( probes.probes ):
            
        t_0 = time.time()

        # -- check if the probe is at the ridge or valley
        xyz = probe.xyz
        if abs(xyz[1]) < 0.001 and abs(xyz[2]) < 0.001: os.chdir( dirs.pp_psd_ridge ) ; n_ridge += 1
        elif abs(xyz[1] + h) < 0.001 :                  os.chdir( dirs.pp_psd_valley ); n_valley += 1
        else:                                           os.chdir( dirs.pp_psd_others ); n_others += 1
        
        probedata = ProbeData( prb_data[i], withT=prb_withT )
        probedata.cleandata( t_start=20.0 )
        probedata.get_fluc( 'p' )
        probedata.pre_multi_psd( 'p_fluc', n_seg=8, overlap=0.5 )
        probedata.write_psd()
        
        t_1 = time.time()
        print(col.fg.green,f"{float(i)/n_data*100:6.2f}%",col.reset,end='')
        print(f" PSD of probe {prb_data[i][-9:-4]} took {t_1-t_0:8.2f}s.")
    
    print(f"\n ridge:{n_ridge}, valley:{n_valley}, others:{n_others}.\n")

