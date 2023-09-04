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

from   vista.psd         import ProbeData

from   vista.timer       import timer

from   vista.tools       import get_filelist

from   vista.colors      import colors as col


probe_dir = os.getcwd()

probe_x_dir = probe_dir + '/probe_x'

psd_x_dir = probe_dir + 'psd_x'

filelist = get_filelist( probe_x_dir )

n_probe = len( filelist )


print(f"We are in directory:{probe_dir}\n")

print(f"We have got {n_probe:5d} probes data.\n")


with timer('Finish computing PSD '):
    
    # enter the psd_x_dir
    
    if os.path.exists( psd_x_dir ): os.chdir( psd_x_dir )
    
    else: os.mkdir( psd_x_dir ); os.chdir( psd_x_dir )
    
    # start 
    
    for i, probefile in enumerate( filelist ):
            
            t_0 = time.time()
            
            probe = ProbeData( probefile, withT=False )
            
            probe.cleandata( starttime=20 )
            
            probe.psd( 8, 0.5 )

            probe.write_psd()
            
            t_1 = time.time()
            print(col.fg.green,f"{float(i)/n_probe*100:6.2f}%",col.reset,end='')
            print(f" PSD of probe {probefile[-9:-4]} took {t_1-t_0:8.2f}s.")

