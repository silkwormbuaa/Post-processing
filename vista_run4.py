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

folderpath = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/'

probefile1 = 'probe_00001.dat'

os.chdir(folderpath)

with timer("reading one probe"):
    
    probe1 = ProbeData(probefile1)

#    print(probe1.df)
    
with timer("clean data"):
    
    probe1.cleandata(starttime=20,Nseg=8)

    print("=============")
    
with timer('psd'):
    
    probe1.psd( 8, 0.5 )
    
    print(probe1.nperseg)
    
    print(probe1.noverlap)
    
    print(len(probe1.freq))
    
    print(len(probe1.pprime_psd))
    
    os.chdir(os.pardir)
    with open("psd.dat",'w') as f:
        for i in range(len(probe1.freq)):
            f.write(str('{:<17.8e}'.format(probe1.freq[i])))
            f.write(str('{:<17.8e}'.format(probe1.pprime_psd[i])))
            f.write(str('{:<17.8e}'.format(probe1.St[i])))
            f.write(str('{:<17.8e}'.format(probe1.nd_pprime_fwpsd[i])))
            f.write('\n')

with timer("plot psd"):
    
    probe1.plot_psd()
    