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
                f.write(str('{:<17.8e}'.format(group.xctr)) + ' ')
                f.write(str('{:<17.8e}'.format(group.yctr)) + ' ')
                f.write(str('{:<17.8e}'.format(group.xmin)) + ' ')
                f.write(str('{:<17.8e}'.format(group.ymin)) + ' ')
                f.write(str('{:<17.8e}'.format(group.xmax)) + ' ')
                f.write(str('{:<17.8e}'.format(group.ymax)) + ' ')
                for name in group.zonelist:
                    f.write(str(name) + ' ')
                f.write('\n')    
    return zonegrp 

#%% Read in the zone group info from txt file.
def ReadZonegrp(FoldPath,filename):
    os.chdir(FoldPath)
    os.chdir(os.pardir)
    zonegrp = list()
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
    return zonegrp

#%% Read the blocks in required region
def ReadBlock(zonegrp,FoldPath):
    RegionRange = [-71.75, 0.0, -60.6875, 86.0] #22.8044
#---Get the zones overlaps with the select region
    zone_overlap_indx = []
    layers_yctr = []
    for i in range(len(zonegrp)):
        notOverlap = (zonegrp[i].xmax <= RegionRange[0]) or \
                     (zonegrp[i].xmin >= RegionRange[2]) or \
                     (zonegrp[i].ymax <= RegionRange[1]) or \
                     (zonegrp[i].ymin >= RegionRange[3])
        Overlap = not notOverlap
        if Overlap: 
            zone_overlap_indx.append(i)
            layers_yctr.append(zonegrp[i].yctr)
            print(i,zonegrp[i].zonelist)
#---Group the zones by their y location
    layer_n = len(np.unique(layers_yctr))
    zonels_y = []
    for i in range(layer_n):
        zonels_y.append([])
    for i in range(len(zone_overlap_indx)):
        if (0.0        < zonegrp[zone_overlap_indx[i]].yctr < 1.73695342):
            zonels_y[0].append(zone_overlap_indx[i])
        elif (1.73695342 < zonegrp[zone_overlap_indx[i]].yctr < 5.01031265):
            zonels_y[1].append(zone_overlap_indx[i])        
        elif (5.01031265 < zonegrp[zone_overlap_indx[i]].yctr < 11.1790910):
            zonels_y[2].append(zone_overlap_indx[i])
        elif (11.1790910 < zonegrp[zone_overlap_indx[i]].yctr < 22.8044042):
            zonels_y[3].append(zone_overlap_indx[i])   
        elif (22.8044042 < zonegrp[zone_overlap_indx[i]].yctr < 44.7127788):
            zonels_y[4].append(zone_overlap_indx[i])
        elif (44.7127788 < zonegrp[zone_overlap_indx[i]].yctr < 86.0000000):
            zonels_y[5].append(zone_overlap_indx[i])   
    print(zonels_y)
#    zonels_y   = [[76,77,78],[134,135],[163],[179],[195],[211]]
         
#---Read in the szplt file with tecplot
    FileList = sorted(GetFileList(FoldPath))
    dataset  = tp.data.load_tecplot_szl(FileList)
    zone     = dataset.zone

#---Manipulate the block data, read in variables, reduce dimension
#    for i in range(len(zonels_y)):
    for i in range(len(zone_overlap_indx)):
        zone_span_n = len(zonegrp[zone_overlap_indx[i]].zonelist)
        for j in range(zone_span_n):
            zonename = zonegrp[zone_overlap_indx[i]].zonelist[j]
            if j==0:
                Nxi, Nyi, Nzi = dataset.zone(zonename).dimensions
                x      = zone(zonename).values('x').as_numpy_array() 
                y      = zone(zonename).values('y').as_numpy_array() 
                u      = zone(zonename).values('<u>').as_numpy_array()
                ugrp   = u
            else:
                u      = zone(zonename).values('<u>').as_numpy_array()
                
                ugrp   = np.vstack((ugrp,u))
#---Span-wise average variables for zone_group i
        if np.ndim(ugrp) == 1:
            umean      = ugrp
        else:
            umean      = np.mean(ugrp,  axis = 0)
#---Reshape the list of variables into matrix
        x = np.unique(x)
        y = np.unique(y)
        umean   = np.reshape(umean,  (Nyi,Nxi))
