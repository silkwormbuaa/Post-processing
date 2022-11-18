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

from   vista_PSD         import * 

from   timer             import timer

from   vista_tools       import *


folderpath = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/probe_x'

outpath = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/psd_x'

os.chdir(folderpath)

filelist = get_filelist( folderpath )

n_probe = len(filelist)

print( "we have got %5d probes data"%( n_probe ) )

os.chdir(outpath)

#with timer('PSD for one case'):
#    
#    for i, probefile in enumerate(filelist):
#        
#        with timer('%5.3f%% Get psd for %s probe'%(float(i)/n_probe*100,probefile[-9:-4])):
#            
#            probe = ProbeData( probefile )
#            
#            probe.cleandata( starttime=20 )
#            
#            probe.psd( 8, 0.5 )
#
#            probe.write_psd()
        

with timer("reading one probe"):
    
    probe1 = ProbeData( filelist[0] )
    
with timer("clean data"):
    
    probe1.cleandata(starttime=20)
    
with timer('psd'):
    
    probe1.psd( 8, 0.5 )
    
    print(probe1.nperseg)
    
    print(probe1.df['pprime'])
    
    os.chdir(os.pardir)
    
#    probe1.write_psd()

