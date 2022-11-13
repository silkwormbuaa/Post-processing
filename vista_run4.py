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
    
    print(probe1.df)
    
    print(probe1.frequency)
    
    print(probe1.meanp)

