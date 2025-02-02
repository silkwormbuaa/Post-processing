#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   probe_250120.py
@Time    :   2025/02/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.probe      import ProbeFile, Probe


def WriteProbe():

    outfile  = '/home/wencanwu/my_simulation/STBLI_low_Re/250120/inca_probes.inp'
    probes = ProbeFile()
#    probes.read( fname )
#    probes.show()

# --- probes at divergent line

    xs = np.arange( -53.6, 115, 0.43382353 )
    y  = 0.00001
    z  = 0.00001
    
    for x in xs:
        
        probe = Probe()
        probe.assign( 'POINT', 5, [x,y,z], 'ALL' )
        probes.probes.append( probe ) 
    
    print(f"Number of probes at the ridges: {len(probes.probes)}. [0,{len(probes.probes)-1}]")
    n = len(probes.probes)
    
# --- probes at the convergent line

    z  = 2.60001
    y  = 0.00001

    for x in xs:
        
        probe = Probe()
        probe.assign( 'POINT', 5, [x,y,z], 'ALL' )
        probes.probes.append( probe ) 
    
    print(f"Number of probes at the valleys: {len(probes.probes)-n}. [{n},{len(probes.probes)-1}]")
    n = len(probes.probes)

    
# --- a vertical line at (x-x_imp)/delta = -5.0

    x = 24.4
    z = 0.001
    ys = np.arange( 0.1, 15.6, 0.2)
    
    for y in ys:
        probe = Probe()
        probe.assign( 'POINT', 5, [x,y,z], 'ALL' )
        probes.probes.append( probe )

    print(f"Number of probes at the vertical line: {len(probes.probes)-n}. [{n},{len(probes.probes)-1}]")
    
    probes.write( outfile )


if __name__ == "__main__":

    WriteProbe()