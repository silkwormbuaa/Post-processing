# -*- coding: utf-8 -*-
'''
@File    :   ReadIn.py
@Time    :   2022/07/06 16:25:13
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Class and functions related to read in data
'''
import tecplot as tp
import os
import numpy as np
from os import path
from timer import timer
#%% #  Define a class ZoneGroup, this class stores the names of zones 
    #  sharing the same location.
class ZoneGroup:
    def __init__(self,xctr,yctr):
        self.xctr = xctr
        self.yctr = yctr
        self.zonelist = []
    xmin     = 0.0
    ymin     = 0.0
    xmax     = 0.0
    ymax     = 0.0
#%%  Get all files list within a directory(contains sub-folders)    
def GetFileList(FoldPath):
    FileList = []
    for home, dirs, files in os.walk(FoldPath):
        for filename in files:
            # Filelist contains the whole path
            FileList.append(os.path.join(home,filename))
    return FileList
#%%  get the zone groups(same x-y location) from statistic data.
def GetZonegrp(FoldPath):
    with timer("Read in statistic data"):
        # get the all szplt file's path+name            
        FileList = sorted(GetFileList(FoldPath))
        print(np.size(FileList))
        # read in szplt file
        dataset  = tp.data.load_tecplot_szl(FileList)
        # VarList  = [v.name for v in dataset.variables()]
        zone     = dataset.zone
        # zonegrp is a list of instances of class ZoneGroup        
        zonegrp  = list()
        # zonectr is a list contains the location of zone centers              
        zonectr  = []
        for i in range(np.size(FileList)):
            zonename = zone(i).name
            # find the range of a zone
            xmin = round(dataset.zone(zonename).values('x').min(),7)
            ymin = round(dataset.zone(zonename).values('y').min(),7)
            xmax = round(dataset.zone(zonename).values('x').max(),7)
            ymax = round(dataset.zone(zonename).values('y').max(),7)
            xctr = 0.5 * (xmin + xmax)
            yctr = 0.5 * (ymin + ymax)
            # if the zone's center is in the list, add this zone's 
            # name to the already existed zone group.
            if [xctr,yctr] in zonectr:
                indx = zonectr.index([xctr,yctr])
                zonegrp[indx].zonelist.append(zonename)
            # if the zone is not in the list, add a new zone group 
            # to the zone group list, then add this zone to this
            # new zone group    
            else:
                zonegrp.append(ZoneGroup(xctr,yctr))            
                zonectr.append([xctr,yctr])
                indx = zonectr.index([xctr,yctr])
                zonegrp[indx].zonelist.append(zonename)
                zonegrp[indx].xmin = xmin
                zonegrp[indx].ymin = ymin
                zonegrp[indx].xmax = xmax
                zonegrp[indx].ymax = ymax
        # sort the zone center coordinates firstly based on y, then x        
        zonectr.sort(key=lambda x: x[1]) 
        # sort the zonegrp by location ang output the groups
        zonegrp.sort(key=lambda x: (x.yctr,x.xctr))     
        # print the zone groups info
        os.chdir(FoldPath)
        os.chdir(os.pardir)
        with open("zonelist.dat","w") as f:
            for group in zonegrp:
                print(group.zonelist)
                f.write(str('{:<16.6e}'.format(group.xctr)) + ' ')
                f.write(str('{:<16.6e}'.format(group.yctr)) + ' ')
                f.write(str('{:<16.6e}'.format(group.xmin)) + ' ')
                f.write(str('{:<16.6e}'.format(group.ymin)) + ' ')
                f.write(str('{:<16.6e}'.format(group.xmax)) + ' ')
                f.write(str('{:<16.6e}'.format(group.xmax)) + ' ')
                for name in group.zonelist:
                    f.write(str(name) + ' ')
                f.write('\n')    
    return zonegrp 

#%% Read in the zone group info from txt file.
def ReadZonegrp(FoldPath,filename):
    os.chdir(FoldPath)
    os.chdir(os.pardir)
    zonegrp = list()
    with timer("Read in zone(block) group info from txt file"):
        with open(filename) as f:
            i = 0
            for line in f.readlines():
                cleanline = line.strip().split()
                xctr = float(cleanline[0])
                yctr = float(cleanline[1])
                zonegrp.append(ZoneGroup(xctr,yctr))
                zonegrp[i].xmin = float(cleanline[2])
                zonegrp[i].ymin = float(cleanline[3])
                zonegrp[i].xmax = float(cleanline[4])
                zonegrp[i].ymax = float(cleanline[5])
                for j in range(6,14):
                    zonegrp[i].zonelist.append(cleanline[j])
                i = i + 1

                