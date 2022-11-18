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

#probefile = 'probe_00142.dat'
#Lsep      = 13.12627403

probefile = 'probe_00135.dat'
Lsep      = 11.340638

with timer("reading one probe"):
    
    probe1 = ProbeData( probefile )
    
with timer("clean data"):
    
    probe1.cleandata(starttime=20)
    
with timer('psd'):
    
    os.chdir(os.pardir)
    
    probe1.psd( 8, 0.5 )
    
    probe1.plot_psd( Lsep )
    
#    print(probe1.nperseg)
    
#    print(probe1.df['pprime'])
    
    
    
#    probe1.write_psd()

