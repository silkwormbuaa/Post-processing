#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   probe_241030.py
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
from   vista.tools      import define_wall_shape


def WriteProbe():

    outfile  = '/home/wencanwu/my_simulation/STBLI_mid_Re/241030/test.inp'
    casecode = '241030'
    y_valley = -0.13679
    probes = ProbeFile()
#    probes.read( fname )
#    probes.show()

# --- probes at the ridges

    xs = np.arange( -53.6, 115, 0.43382353 )
    y  = 0.00001
    z  = 0.00001
    
    for x in xs:
        
        probe = Probe()
        probe.assign( 'POINT', 20, [x,y,z], 'ALL' )
        probes.probes.append( probe ) 
    
    print(f"Number of probes at the ridges: {len(probes.probes)}. [0,{len(probes.probes)-1}]")
    n = len(probes.probes)
    
# --- probes at the valleys

    z  = 0.1625
    y  = y_valley

    for x in xs:
        
        probe = Probe()
        probe.assign( 'POINT', 20, [x,y,z], 'ALL' )
        probes.probes.append( probe ) 
    
    print(f"Number of probes at the valleys: {len(probes.probes)-n}. [{n},{len(probes.probes)-1}]")
    n = len(probes.probes)
    
# --- probes before the interaction region

    x = -53.6
    zs = np.arange( -5.2, 0.065, 0.065)
    ys = define_wall_shape( zs, casecode=casecode, write=False )
    
    for z,y in zip(zs,ys):
        
        probe = Probe()
        probe.assign( 'POINT', 20, [x,y,z], 'ALL' )
        probes.probes.append( probe )
    
    print(f"Number of probes before the interaction region: {len(probes.probes)-n}. [{n},{len(probes.probes)-1}]")
    n = len(probes.probes)
    
# --- probes after the interaction region

    x = 76.4
    
    for z,y in zip(zs,ys):
        
        probe = Probe()
        probe.assign( 'POINT', 20, [x,y,z], 'ALL' )
        probes.probes.append( probe )

    print(f"Number of probes after the interaction region: {len(probes.probes)-n}. [{n},{len(probes.probes)-1}]")
    n = len(probes.probes)
    
# --- a vertical line at (x-x_imp)/delta = -5.0

    x = 24.4
    z = 0.001
    ys = np.arange( 0.1, 15.6, 0.2)
    
    for y in ys:
        probe = Probe()
        probe.assign( 'POINT', 20, [x,y,z], 'ALL' )
        probes.probes.append( probe )

    print(f"Number of probes at the vertical line: {len(probes.probes)-n}. [{n},{len(probes.probes)-1}]")
    
    probes.write( outfile )


if __name__ == "__main__":

    WriteProbe()