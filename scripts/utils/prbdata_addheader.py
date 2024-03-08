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

source_dir  = '/media/wencanwu/Seagate Expansion Drive1/temp/240211/probes' 
target_dir  = '/media/wencanwu/Seagate Expansion Drive1/temp/240210/probes'
outputpath  = '/media/wencanwu/Seagate Expansion Drive1/temp/test/output'

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
