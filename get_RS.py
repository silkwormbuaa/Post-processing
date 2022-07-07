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
    uulst= []
    vvlst= []
    wwlst= []
    uvlst= []
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
                uu = zone(zonename).values('<u`u`>').as_numpy_array()
                uv = zone(zonename).values('<u`v`>').as_numpy_array()
                vv = zone(zonename).values('<v`v`>').as_numpy_array()
                ww = zone(zonename).values('<w`w`>').as_numpy_array()
                rho = zone(zonename).values('<rho>').as_numpy_array()
                mu = zone(zonename).values('<mu>').as_numpy_array()
                ugrp   = u
                uugrp  = uu
                uvgrp  = uv
                vvgrp  = vv
                wwgrp  = ww
                rhogrp = rho                
                mugrp  = mu
            else:
                u   = zone(zonename).values('<u>').as_numpy_array()
                uu = zone(zonename).values('<u`u`>').as_numpy_array()
                uv = zone(zonename).values('<u`v`>').as_numpy_array()
                vv = zone(zonename).values('<v`v`>').as_numpy_array()
                ww = zone(zonename).values('<w`w`>').as_numpy_array()                
                rho = zone(zonename).values('<rho>').as_numpy_array()
                mu  = zone(zonename).values('<mu>').as_numpy_array()
                # ugrp stores u for all spanwise zones(u of the group)
                ugrp = np.vstack((ugrp,u))
                uugrp= np.vstack((uugrp,uu))
                uvgrp= np.vstack((uvgrp,uv))
                vvgrp= np.vstack((vvgrp,vv))
                wwgrp= np.vstack((wwgrp,ww))
                rhogrp = np.vstack((rhogrp,rho))
                mugrp  = np.vstack((mugrp,mu))
        # spanwise average the u for the zone group
        # when ugrp has only one dimension, no need to to average
        if np.ndim(ugrp) == 1:
            umean = ugrp
            uumean  = uugrp
            uvmean  = uvgrp
            vvmean  = vvgrp 
            wwmean  = wwgrp
            rhomean = rhogrp
            mumean  = mugrp
        else:                     
            umean = np.mean(ugrp,axis = 0)
            uumean = np.mean(uugrp,axis = 0)
            uvmean = np.mean(uvgrp,axis = 0)
            vvmean = np.mean(vvgrp,axis = 0)
            wwmean = np.mean(wwgrp,axis = 0)
            rhomean = np.mean(rhogrp,axis = 0)
            mumean  = np.mean(mugrp,axis = 0)
        x = np.unique(x)
        y = np.unique(y)
        # notice the order of u
        umean = np.reshape(umean,(Nyi,Nxi))
        uumean  = np.reshape(uumean,(Nyi,Nxi))
        uvmean  = np.reshape(uvmean,(Nyi,Nxi))
        vvmean  = np.reshape(vvmean,(Nyi,Nxi))
        wwmean  = np.reshape(wwmean,(Nyi,Nxi))
        rhomean = np.reshape(rhomean,(Nyi,Nxi))
        mumean  = np.reshape(mumean,(Nyi,Nxi))
        
        # search for the points for interpolation
        for k in range(np.size(x)):
            if x[k] <= line_loc[0] < x[k+1]:
                x_i = k 
                break
        # do interpolation
        ## "if" here to avoid picking up the point at zone interface twice
        if i < (len(zone_intersect_i)-1):
            y_i_max = np.size(y)-1
        else:
            y_i_max = np.size(y)             
        for y_i in range(y_i_max):
            local_u = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                                umean[y_i][x_i],umean[y_i][x_i+1])
            local_uu  = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                                uumean[y_i][x_i],uumean[y_i][x_i+1])   
            local_uv  = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                                uvmean[y_i][x_i],uvmean[y_i][x_i+1])
            local_vv  = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                                vvmean[y_i][x_i],vvmean[y_i][x_i+1])
            local_ww  = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                                wwmean[y_i][x_i],wwmean[y_i][x_i+1])                                                   
            local_rho = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                                rhomean[y_i][x_i],rhomean[y_i][x_i+1])
            local_mu  = GetInterpo(line_loc[0],x[x_i],x[x_i+1],\
                                mumean[y_i][x_i],mumean[y_i][x_i+1])    
            # if the point is in the range of line(y direction)        
            if (y[y_i] >= line_loc[1]) & (y[y_i] <= line_loc[3]) :
                xlst.append(line_loc[0])
                #xlst.append(x[x_i])
                ylst.append(y[y_i])            
                ulst.append(local_u)
                uulst.append(local_uu)
                uvlst.append(local_uv)
                vvlst.append(local_vv)
                wwlst.append(local_ww)
                rholst.append(local_rho)
                mulst.append(local_mu)
    # get the friction velocity and viscous length scale at the wall
    dudy = abs((ulst[1]-ulst[0])/(ylst[1]-ylst[0]))        
    tau  = mulst[0] * dudy
    u_tau= np.sqrt(tau/rholst[0])
    lnu  = mulst[0]/rholst[0]/u_tau                 
    # plot curve
    os.chdir(OutFile)
    out_file = "RS_x_" + str(line_loc[0]) + "flat.dat"
    with open(out_file,"w") as f:
        for i in range(len(ylst)):
            yplus  = (ylst[i] - ylst[0])/lnu
            uunor = np.sqrt(abs(uulst[i]))*np.sqrt(rholst[i]/rholst[0])/u_tau
            vvnor = np.sqrt(vvlst[i])*np.sqrt(rholst[i]/rholst[0])/u_tau
            wwnor = np.sqrt(wwlst[i])*np.sqrt(rholst[i]/rholst[0])/u_tau
            uvnor = np.sqrt(abs(uvlst[i])*rholst[i]/rholst[0]/u_tau**2)*uvlst[i]/abs(uvlst[i])
            f.write(str(yplus)+' '+str(uunor)+' '+str(uvnor)+' '+ \
                    str(vvnor)+' '+str(wwnor)+'\n')
    
    
#    plt.show()
        


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

line_loc = [-58.25, 0, -58.25, 30.0]    
GetLine(line_loc,zonegrp)
#line_loc = [-57.0, 0.0, -57.0, 30.0]    
#line_loc = [-97.25, -0.52, -97.25, 30.0]
#GetLine(line_loc,zonegrp)
#line_loc = [-80.0, 0.0, -80.0, 10.0]    
#GetLine(line_loc,zonegrp)
#line_loc = [-60.0, 0.0, -60.0, 10.0]    
#GetLine(line_loc,zonegrp)
#line_loc = [-40.0, 0.0, -40.0, 10.0]    
#GetLine(line_loc,zonegrp)
#line_loc = [-20.0, 0.0, -20.0, 10.0]    
#GetLine(line_loc,zonegrp)
#line_loc = [0.0, 0.0, 0.0, 10.0]    
#GetLine(line_loc,zonegrp)


#        print(dataset.num_solution_times)
#        print(dataset.VariablesNamedTuple)
#        print(dataset.num_variables)
#        print(dataset.num_zones)
#        print(dataset.title)
#        print(dataset.variable_names)
