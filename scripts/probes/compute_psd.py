#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   compute_psd.py
@Time    :   2022/11/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Scripts for getting PSD
'''

import os
import sys
import time
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.probe       import ProbeFile
from   vista.timer       import timer
from   vista.params      import Params
from   vista.colors      import colors as col
from   vista.directories import Directories
from   vista.material    import get_visc
from   vista.directories import create_folder
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

dirs   = Directories( os.getcwd() )
probes = ProbeFile( dirs.set_prb )

# -- enter working directory and read relevant parameters 

os.chdir( dirs.prb_dir )

prb_data = dirs.probes
n_data = len( prb_data )
print(f"We are in directory:{dirs.prb_dir}\n")
print(f"We have got {n_data:5d} probes data.\n")

params     = Params( dirs.case_para_file )
h          = params.H
u_ref      = params.u_ref
rho_ref    = params.rho_ref
delta      = params.delta_0
Re_ref     = params.Re_ref
visc_law   = params.visc_law
prb_withT  = params.prb_withT
p_dyn      = 0.5 * rho_ref * u_ref**2

# -- check if number of probe data and number of probes are consistent

if n_data != len( probes.probes ):
    raise ValueError("Number of probe data and number of probes are not consistent.")

# -- check if probe location and classify them to either ridge/valley or others

create_folder( dirs.pp_psd_ridge  )
create_folder( dirs.pp_psd_valley )
create_folder( dirs.pp_psd_others )

# -- start computing PSD 

with timer('Computing PSD '):
    
    # start 
    
    n_ridge = 0 ; n_valley = 0 ; n_others = 0
    
    for i in range( 0, len(probes.probes) ):
            
        t_0 = time.time()
        
        probe = probes.probes[i]

        # -- check if the probe is at the ridge or valley
        xyz = probe.xyz
        if abs(xyz[1]) < 0.001 and abs(xyz[2]) < 0.001: os.chdir( dirs.pp_psd_ridge ) ; n_ridge += 1
        elif abs(xyz[1] + h) < 0.001 :                  os.chdir( dirs.pp_psd_valley ); n_valley += 1
        else:                                           os.chdir( dirs.pp_psd_others ); n_others += 1
        
        probedata = ProbeData( prb_data[i], withT=prb_withT )
        probedata.cleandata( t_start=20.0 )
        
        # -- compute cf first
        
        if os.getcwd() == dirs.pp_psd_ridge:
            walldist = probedata.xyz[1]
        elif os.getcwd() == dirs.pp_psd_valley:
            walldist = abs( h + probedata.xyz[1] )
        else: walldist = probedata.xyz[1]
        
        # in case the old cases that not having T in the probe data
        
        try:
            ts  = np.array( probedata.df['T'] )
        except:
            ts  = np.array(probedata.df['p']/probedata.df['rho']/287.0508571)
        
        u       = np.array( probedata.df['u'] )
        mu      = get_visc( ts, Re_ref, law= visc_law )
        cf      = mu*u/walldist/p_dyn
        
        probedata.add_data( 'cf', cf )
        probedata.get_fluc( ['cf','p'] )
        probedata.compute_psd(   'cf_fluc', n_seg=8, overlap=0.5 )
        probedata.compute_psd(   'p_fluc',  n_seg=8, overlap=0.5 )
        probedata.pre_multi_psd( 'cf_fluc', n_seg=8, overlap=0.5 )
        probedata.pre_multi_psd( 'p_fluc',  n_seg=8, overlap=0.5 )
        probedata.write_psd()
        
        t_1 = time.time()
        print(col.fg.green,f"{float(i)/n_data*100:6.2f}%",col.reset,end='')
        print(f" PSD of probe {prb_data[i][-9:-4]} took {t_1-t_0:8.2f}s.")
    
    print(f"\n ridge:{n_ridge}, valley:{n_valley}, others:{n_others}.\n")

# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
sys.stdout.flush()   
