# -*- coding: utf-8 -*-
'''
@File    :   BLprofile.py
@Time    :   2022/07/06 15:40:10
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Read in statistic data and extract data on a line
             to get profile data.   
'''
from pyexpat.errors import XML_ERROR_INCOMPLETE_PE
import tecplot as tp
import os
from os import path
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
from timer import timer
from ReadIn import *

#%% Get linear interpolation 
def GetInterpo(x,x1,x2,y1,y2):
    slope = (y2-y1)/(x2-x1)
    y = slope * (x-x1) + y1
    return y
#%% Get velocity gradient magnitude at the wall
def GetDudy(u1,u2,y1,y2,opt=2):
    if opt == 1:            #  linear gradient
        dudy = abs(u1/y1)
    if opt == 2:            #  second order accuracy
        temp1 = u1*y2*y2 - u2*y1*y1
        temp2 = y1*y2*(y2-y1)
        dudy = abs(temp1/temp2)
    return dudy 
#%% Get the profile along a normal line
def GetLine(line_loc,zonegrp,FoldPath,OutPath):
    # line_loc = [x1, y1,  x2, y2]
    # set an empty list for storing the intersected zone
    # only vertical line works so far
    ylst = []
    xlst = []
    ulst = []
    zone_intersect_i = []
    for i in range(len(zonegrp)):
        line_intersected =  (line_loc[0] >= zonegrp[i].xmin) and \
                            (line_loc[0] <  zonegrp[i].xmax) and \
                            (line_loc[1] <  zonegrp[i].ymax) and \
                            (line_loc[3] >  zonegrp[i].ymin)                        
        if line_intersected:
            zone_intersect_i.append(i)
    # HAVE SELECTED ZONES FOR AVERAGING| START TO READ IN DATA    
    FileList = sorted(GetFileList(FoldPath))
    dataset  = tp.data.load_tecplot_szl(FileList)
    zone     = dataset.zone
    for i in range(len(zone_intersect_i)):
        # for each zone group, x,y DO NOT need to be averaged
        # other variables need to be averaged
        
        # firstly get the total number of zones in a zone group
        zone_span_n = len(zonegrp[zone_intersect_i[i]].zonelist)
        for j in range(zone_span_n):
            zonename = zonegrp[zone_intersect_i[i]].zonelist[j] 
            # read in x and y for the first time  
            if j==0:
                Nxi, Nyi, Nzi = dataset.zone(zonename).dimensions
                x = zone(zonename).values('x').as_numpy_array() 
                y = zone(zonename).values('y').as_numpy_array() 
                u = zone(zonename).values('<u>').as_numpy_array()
                ugrp = u
            else:
                u = zone(zonename).values('<u>').as_numpy_array()
                # ugrp stores u for all spanwise zones(u of the group)
                ugrp = np.vstack((ugrp,u))
        # spanwise average the u for the zone group
        # when ugrp has only one dimension, no need to to average
        if np.ndim(ugrp) == 1:
            umean = ugrp
        else:                     
            umean = np.mean(ugrp,axis = 0)
        x = np.unique(x)
        y = np.unique(y)
        # notice the order of u
        umean = np.reshape(umean,(Nyi,Nxi))
        
        # search for the points for interpolation
        for k in range(np.size(x)):
            if x[k] <= line_loc[0] < x[k+1]:
                x_i = k
                break
        # do interpolation
        if i < (len(zone_intersect_i)-1):
            y_i_max = np.size(y)-1
        else:
            y_i_max = np.size(y)             
        for y_i in range(y_i_max):
            xlst.append(line_loc[0])
            ylst.append(y[y_i])
            local_u = GetInterpo(line_loc[0],x[k],x[k+1],\
                                umean[y_i][k],umean[y_i][k+1])
            ulst.append(local_u)
    # plot curve
    os.chdir(OutPath)
    out_file = "u_x_" + str(line_loc[0]) + "flat.dat"
    with open(out_file,"w") as f:
        for i in range(len(ylst)):
            f.write(str(ylst[i])+' '+str(ulst[i])+'\n')