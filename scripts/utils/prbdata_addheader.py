#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   prbdata_addheader.py
@Time    :   2024/03/08 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.tools       import get_filelist
from   vista.directories import Directories
from   vista.probe       import ProbeFile
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.colors      import colors    as col


# =============================================================================

# opt = 1        # copy header from source file to target file
opt = 2          # add header from inca_probe.inp file ()  

# =============================================================================

if opt == 1:

    source_dir  = '/media/wencan/Expansion/temp/240211/probes' 
    target_dir  = '/media/wencan/Expansion/temp/240210/probes'
    outputpath  = '/media/wencan/Expansion/temp/test/output'

    sourcefiles = get_filelist( source_dir, 'probe_' )
    targetfiles = get_filelist( target_dir, 'probe_' )

    os.chdir( outputpath )
    for i, sourcefile in enumerate( sourcefiles ):
        
        with open( sourcefile, 'r' ) as f:
            
            lines = f.readlines()
        
        with open( targetfiles[i], 'r' ) as f:
            
            lines_target = f.read()
        
        with open( targetfiles[i].split('/')[-1], 'w') as f:
                
            f.write( lines[0] )
            f.write( lines_target )
            
        print( f"File {i} is done!")


elif opt == 2:
    
    casepath = '/media/wencan/Expansion/temp/231124'
    
    case = Directories( casepath )
    
    probes = ProbeFile( case.set_prb )
    grid   = GridData( case.grid )
    grid.read_grid()
    
    prbdatafiles = get_filelist( case.prb_dir, 'probe_' )
    
    # check probe numbers
    
    if len( prbdatafiles ) != len( probes.probes ):
        raise ValueError("Probe number and data number are not matched!")
    
    # add header to probe data file
    with timer("Add header to probe data file"):
        
        for i in range( 0, len(prbdatafiles) ):
            
            prbdata = prbdatafiles[i]
            
            xyz_prb = probes.probes[i].xyz
            xyz_grd = grid.find_probe_xyz( xyz_prb )
            
            print( f"Probe {i+1} {xyz_prb} has actual grid coordinate {xyz_grd}")
            
            with open( prbdata, 'r+' ) as f:
                
                content = f.read()
                
                f.seek(0)
                f.truncate()
                
                string = f' x = {xyz_grd[0]}, y = {xyz_grd[1]}, z = {xyz_grd[2]}\n'
                f.write( string )
                f.write( content )
                
            print(col.fg.green,f"Progress {100*i/len(prbdatafiles):7.2f}%",col.reset,f"Probe {i+1} is done!")
