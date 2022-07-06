# -*- coding: utf-8 -*-
"""
Original code credit to Weibo
"""
from pyexpat.errors import XML_ERROR_INCOMPLETE_PE
import tecplot as tp
#from tecplot.constant import *
import pandas as pd
import os
from os import path
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
from timer import timer
from time import time
import logging as log
from glob import glob

log.basicConfig(level=log.INFO)

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
#%%   Get all files list within a directory(contains sub-folders)    
def GetFileList(path):
    FileList = []
    for home, dirs, files in os.walk(path):
        for filename in files:
            # Filelist contains the whole path
            FileList.append(os.path.join(home,filename))
    return FileList
#%% Get linear interpolation 
def GetInterpo(x,x1,x2,y1,y2):
    slope = (y2-y1)/(x2-x1)
    y = slope * (x-x1) + y1
    return y
#%% Line intersect with zone and get data along this line(zones)
def GetLine(line_loc,zonegrp):
    # line_loc = [x1, y1,  x2, y2]
    # set an empty list for storing the intersected zone
    # only vertical line works so far
    ylst = []
    xlst = []
    ulst = []
    Tlst = []
    rholst = []
    mulst  = []    
    zone_intersect_i = []
    for i in range(len(zonegrp)):
        line_intersected =  (line_loc[0] >= zonegrp[i].xmin) and \
                            (line_loc[0] <  zonegrp[i].xmax) and \
                            (line_loc[1] <  zonegrp[i].ymax) and \
                            (line_loc[3] >  zonegrp[i].ymin)                        
        if line_intersected:
            zone_intersect_i.append(i)
    # HAVE SELECTED ZONES FOR AVERAGING| START TO READ IN DATA    

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
                rho = zone(zonename).values('<rho>').as_numpy_array()
                mu  = zone(zonename).values('<mu>').as_numpy_array()
                T   = zone(zonename).values('<T>').as_numpy_array()
                ugrp = u
                rhogrp = rho                
                mugrp  = mu 
                Tgrp   = T               
            else:
                u = zone(zonename).values('<u>').as_numpy_array()
                rho = zone(zonename).values('<rho>').as_numpy_array()
                mu  = zone(zonename).values('<mu>').as_numpy_array()
                T   = zone(zonename).values('<T>').as_numpy_array()                
                # ugrp stores u for all spanwise zones(u of the group)
                ugrp = np.vstack((ugrp,u))
                rhogrp = np.vstack((rhogrp,rho))
                mugrp  = np.vstack((mugrp,mu))
                Tgrp   = np.vstack((Tgrp,T))
        # spanwise average the u for the zone group
        # when ugrp has only one dimension, no need to to average
        if np.ndim(ugrp) == 1:
            umean = ugrp
            rhomean = rhogrp
            mumean  = mugrp
            Tmean   = Tgrp            
        else:                     
            umean = np.mean(ugrp,axis = 0)
            rhomean = np.mean(rhogrp,axis = 0)
            mumean  = np.mean(mugrp,axis = 0)  
            Tmean  = np.mean(Tgrp,axis = 0)          
        x = np.unique(x)
        y = np.unique(y)
        # notice the order of u
        umean = np.reshape(umean,(Nyi,Nxi))
        rhomean = np.reshape(rhomean,(Nyi,Nxi))
        mumean  = np.reshape(mumean,(Nyi,Nxi))  
        Tmean  = np.reshape(Tmean,(Nyi,Nxi))   
           
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

            local_u = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                                umean[y_i][x_i],umean[y_i][x_i+1])
            local_rho = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                            rhomean[y_i][x_i],rhomean[y_i][x_i+1])
            local_mu  = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                            mumean[y_i][x_i],mumean[y_i][x_i+1])
            local_T   = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                            Tmean[y_i][x_i],Tmean[y_i][x_i+1])    
            if (y[y_i] >= line_loc[1]) & (y[y_i] <= line_loc[3]) :   
                xlst.append(line_loc[0])
                ylst.append(y[y_i])                                             
                ulst.append(local_u)
                rholst.append(local_rho)
                mulst.append(local_mu)
                Tlst.append(local_T)  
    # get the friction velocity and viscous length scale at the wall
    dudy = abs((ulst[1]-ulst[0])/(ylst[1]-ylst[0]))        
    tau  = mulst[0] * dudy
    u_tau= np.sqrt(tau/rholst[0])
    lnu  = mulst[0]/rholst[0]/u_tau                            
    # plot curve
    os.chdir(OutFile)
    out_file = "T_x_" + str(line_loc[0]) + "flat.dat"
    with open(out_file,"w") as f:
        for i in range(len(ylst)):
            yplus  = (ylst[i] - ylst[0])/lnu
            yplus2 = ylst[i]/lnu            
            f.write(str(yplus)+' '+str(Tlst[i])+' '+str(yplus2)+' '+str(ylst[i])+'\n')
    
    
        


#%% Read plt data from INCA
FoldPath = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/TP_stat"
OutFile  = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/DataPost/"
#FoldPath = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/TP_stat"
#OutFile  = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/DataPost/"
'''
VarList  = ['x',                'y',                    '<u>',                  \
            '<v>',              '<w>',                  '<rho>',                \
            '<rho*E>',          '<p>',                  '<T>',                  \
            '<mu>',             '||grad(<velocity>)||', '||curl(<velocity>)||', \
            'div(<velocity>)',  '||shear(<velocity>)||','||grad(<T>)||',        \
            '<u`u`>',           '<u`v`>',               '<u`w`>',               \
            '<u`rho`>',         '<u`rho*E`>',           '<u`p`>',               \
            '<u`T`>',           '<u`mu`>',              '<v`v`>',               \
            '<v`w`>',           '<v`rho`>',             '<v`rho*E`>',           \
            '<v`p`>',           '<v`T`>',               '<v`mu`>',              \
            '<w`w`>',           '<w`rho`>',             '<w`rho*E`>',           \
            '<w`p`>',           '<w`T`>',               '<w`mu`>',              \
            '<rho`rho`>',       '<rho`rho*E`>',         '<rho`p`>',             \
            '<rho`T`>',         '<rho`mu`>',            '<rho*E`rho*E`>',       \
            '<rho*E`p`>',       '<rho*E`T`>',           '<rho*E`mu`>',          \
            '<p`p`>',           '<p`T`>',               '<p`mu`>',              \
            '<T`T`>',           '<T`mu`>',              '<mu`mu`>',             \
            'walldist',         'shocksens',            'entropy']
'''
with timer("Read Meanflow data"):
    FileList = sorted(GetFileList(FoldPath))
    print(np.size(FileList))

#   read in szplt file
    dataset  = tp.data.load_tecplot_szl(FileList)
#   print the zone name.
    VarList  = [v.name for v in dataset.variables()]
    zone     = dataset.zone
    zonegrp  = list()
    zonectr  = []
#    for i in range(50):    
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
        # to the zone group list, then add this zone to the new
        # zone group    
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
    for group in zonegrp:
        print(group.zonelist)

line_loc = [-58.25, 0.0, -58.25, 30.0]    
GetLine(line_loc,zonegrp)
line_loc = [-57.00, 0.0, -57.0, 30.0]    
GetLine(line_loc,zonegrp)
line_loc = [-97.25, 0.0, -97.25, 30.0]
GetLine(line_loc,zonegrp)


#        print(dataset.num_solution_times)
#        print(dataset.VariablesNamedTuple)
#        print(dataset.num_variables)
#        print(dataset.num_zones)
#        print(dataset.title)
#        print(dataset.variable_names)
