# -*- coding: utf-8 -*-
'''
@File    :   block.py
@Time    :   2022/09/29 16:49:01
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Making 8 layers of statistics data to one file
'''

import tecplot as tp
import os
import numpy as np
from timer import timer
from ReadIn import * 

#%% Read the szplt files and manipulate data.
def ReadInZone(FoldPath):
    with timer("Read in statistic data"):
        # get the all szplt file's path+name            
        FileList = sorted(GetFileList(FoldPath))
        print(np.size(FileList))
        # read in szplt file
        dataset  = tp.data.load_tecplot_szl(FileList)
    return dataset.zone


FoldPath = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/TP_stat"
#zone = ReadInZone(FoldPath)
#print(type(zone))
FileList = sorted(GetFileList(FoldPath))
dataset  = tp.data.load_tecplot_szl(FileList)
zonegrp = ReadZonegrp(FoldPath,'zonelist.dat')
zone    = dataset.zone
with timer("adding z values to zones"):
    for grp in zonegrp:
        for i in range(len(grp.zonelist)):
            tp.data.operate.execute_equation('{z} = '+str(i*1.3), zones = [zone(grp.zonelist[i])])
print(zone(0).values('z'))


print(os.getcwd())
tp.data.save_tecplot_plt("datawithz.plt")