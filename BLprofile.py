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
import os

import numpy             as np

import tecplot           as tp

from   vista.tools       import get_filelist

#%% Get linear interpolation 
def GetInterpo(x,x1,x2,y1,y2):
    slope = (y2-y1)/(x2-x1)
    y = slope * (x-x1) + y1
    return y

#%% Get velocity gradient magnitude at the wall
def GetDudy(u1,u2,y1,y2,opt=1):
    if opt == 1:            #  linear gradient
        dudy = abs(u1/y1)
    if opt == 2:            #  second order accuracy
        temp1 = u1*y2*y2 - u2*y1*y1
        temp2 = y1*y2*(y2-y1)
        dudy = abs(temp1/temp2)
    return dudy

#%% Get the profile along a normal line
def GetLine(line_loc,zonegrp,FoldPath,OutPath,wall_opt=1):
    # line_loc = [x1, y1,  x2, y2]
    # wall_opt for choosing wall condition, 1="adiabatic",2="isothermal"
    # zonegrp is a list of instances of Class ZoneGroup
    # set an empty list for storing the intersected zone
    # only vertical line works so far
    xlst   = [line_loc[0]]
    ylst   = [line_loc[1]]
    ulst   = [0.0]
    uulst  = [0.0]
    vvlst  = [0.0]
    wwlst  = [0.0]
    uvlst  = [0.0]
    rholst = [0.0]
    mulst  = [0.0]
    Tlst   = [0.0]
    zone_intersect_i = []
    for i in range(len(zonegrp)):
        line_intersected =  (line_loc[0] >= zonegrp[i].xmin) and \
                            (line_loc[0] <  zonegrp[i].xmax) and \
                            (line_loc[1] <  zonegrp[i].ymax) and \
                            (line_loc[3] >  zonegrp[i].ymin)                        
        if line_intersected:
            # store the index of intersected zones
            zone_intersect_i.append(i)  
    # have selected zones intersected | start to read in specific zones   
    FileList = sorted(get_filelist(FoldPath))
    dataset  = tp.data.load_tecplot_szl(FileList)
    zone     = dataset.zone
    for i in range(len(zone_intersect_i)):
        # For each zone group, x,y DO NOT need to be averaged
        # other variables need to be averaged
        # Firstly get the total number of zones in a zone group
        zone_span_n = len(zonegrp[zone_intersect_i[i]].zonelist)
        for j in range(zone_span_n):
            zonename = zonegrp[zone_intersect_i[i]].zonelist[j] 
            # read in x and y for the first time  
            if j==0:
                Nxi, Nyi, Nzi = dataset.zone(zonename).dimensions
                x      = zone(zonename).values('x').as_numpy_array() 
                y      = zone(zonename).values('y').as_numpy_array() 
                u      = zone(zonename).values('<u>').as_numpy_array()
                uu     = zone(zonename).values('<u`u`>').as_numpy_array()
                vv     = zone(zonename).values('<v`v`>').as_numpy_array()
                ww     = zone(zonename).values('<w`w`>').as_numpy_array()
                uv     = zone(zonename).values('<u`v`>').as_numpy_array()                
                rho    = zone(zonename).values('<rho>').as_numpy_array()
                mu     = zone(zonename).values('<mu>').as_numpy_array()
                T      = zone(zonename).values('<T>').as_numpy_array()
                ugrp   = u
                uugrp  = uu
                vvgrp  = vv
                wwgrp  = ww
                uvgrp  = uv
                rhogrp = rho                
                mugrp  = mu
                Tgrp   = T
            else:
                u      = zone(zonename).values('<u>').as_numpy_array()
                uu     = zone(zonename).values('<u`u`>').as_numpy_array()
                vv     = zone(zonename).values('<v`v`>').as_numpy_array()
                ww     = zone(zonename).values('<w`w`>').as_numpy_array()
                uv     = zone(zonename).values('<u`v`>').as_numpy_array()                
                rho    = zone(zonename).values('<rho>').as_numpy_array()
                mu     = zone(zonename).values('<mu>').as_numpy_array()
                T      = zone(zonename).values('<T>').as_numpy_array()                
                # ugrp stores u for all spanwise zones(u of the group)
                ugrp   = np.vstack((ugrp,u))
                uugrp  = np.vstack((uugrp,uu))
                vvgrp  = np.vstack((vvgrp,vv))
                wwgrp  = np.vstack((wwgrp,ww))
                uvgrp  = np.vstack((uvgrp,uv))             
                rhogrp = np.vstack((rhogrp,rho))
                mugrp  = np.vstack((mugrp,mu))
                Tgrp   = np.vstack((Tgrp,T))
        # spanwise average the variables for zone group i
        # when ugrp has only one dimension, no need to to average
        if np.ndim(ugrp) == 1:
            # if here use np.mean, will collapse to a single value
            umean   = ugrp
            uumean  = uugrp
            vvmean  = vvgrp 
            wwmean  = wwgrp
            uvmean  = uvgrp            
            rhomean = rhogrp
            mumean  = mugrp
            Tmean   = Tgrp         
        else:                     
            umean   = np.mean(ugrp,  axis = 0)
            uumean  = np.mean(uugrp, axis = 0)
            vvmean  = np.mean(vvgrp, axis = 0)
            wwmean  = np.mean(wwgrp, axis = 0)
            uvmean  = np.mean(uvgrp, axis = 0)
            rhomean = np.mean(rhogrp,axis = 0)
            mumean  = np.mean(mugrp, axis = 0)
            Tmean   = np.mean(Tgrp,  axis = 0)
        # reshape the list of variables into matrix    
        x = np.unique(x)
        y = np.unique(y)
        # notice the order of u
        umean   = np.reshape(umean,  (Nyi,Nxi))
        uumean  = np.reshape(uumean, (Nyi,Nxi))
        vvmean  = np.reshape(vvmean, (Nyi,Nxi))
        wwmean  = np.reshape(wwmean, (Nyi,Nxi))
        uvmean  = np.reshape(uvmean, (Nyi,Nxi))
        rhomean = np.reshape(rhomean,(Nyi,Nxi))
        mumean  = np.reshape(mumean, (Nyi,Nxi))
        Tmean   = np.reshape(Tmean,  (Nyi,Nxi))
        # search for the points for interpolation
        for k in range(np.size(x)):
            if x[k] <= line_loc[0] < x[k+1]:
                break
        # do interpolation(in x direction)
        # check if it is a top block, if not, yimax = size(y)-1 \
        # to avoid double count overlapped point
        if i < (len(zone_intersect_i)-1):
            y_i_max = np.size(y)-1
        else:
            y_i_max = np.size(y)             
        for y_i in range(y_i_max):
            local_u   = GetInterpo(line_loc[0],    x[k], x[k+1],\
                                   umean[y_i][k],  umean[y_i][k+1])
            local_uu  = GetInterpo(line_loc[0],    x[k], x[k+1],\
                                   uumean[y_i][k], uumean[y_i][k+1])   
            local_vv  = GetInterpo(line_loc[0],    x[k], x[k+1],\
                                   vvmean[y_i][k], vvmean[y_i][k+1])
            local_ww  = GetInterpo(line_loc[0],    x[k], x[k+1],\
                                   wwmean[y_i][k], wwmean[y_i][k+1])
            local_uv  = GetInterpo(line_loc[0],    x[k], x[k+1],\
                                   uvmean[y_i][k], uvmean[y_i][k+1])
            local_rho = GetInterpo(line_loc[0],    x[k], x[k+1],\
                                   rhomean[y_i][k],rhomean[y_i][k+1])
            local_mu  = GetInterpo(line_loc[0],    x[k], x[k+1],\
                                   mumean[y_i][k], mumean[y_i][k+1])
            local_T   = GetInterpo(line_loc[0],    x[k], x[k+1],\
                                   Tmean[y_i][k],  Tmean[y_i][k+1])   
            # if the point is in the range of line(y direction) 
            if (y[y_i] >= line_loc[1]) & (y[y_i] <= line_loc[3]) :
                if y[y_i] == line_loc[1] and wall_opt==2:
                    rholst[0] = local_rho
                    mulst[0]  = local_mu
                    Tlst[0]   = local_T
                    # print("x=",line_loc[0]," ",rholst[0]," ",mulst[0]," ",Tlst[0])
                elif y[y_i] > line_loc[1]:    
                    xlst.append(line_loc[0])
                    #xlst.append(x[x_i])
                    ylst.append(y[y_i])            
                    ulst.append(local_u)
                    uulst.append(local_uu)
                    vvlst.append(local_vv)
                    wwlst.append(local_ww)
                    uvlst.append(local_uv)
                    rholst.append(local_rho)
                    mulst.append(local_mu)
                    Tlst.append(local_T)
    # get the friction velocity and viscous length scale at the wall
    if wall_opt==1:
        rholst[0] = rholst[1]
        mulst[0]  = mulst[1]
        Tlst[0]   = Tlst[1]
    # dudy opt==1, 1st order; opt==2, 2nd order    
    dudy      = GetDudy(ulst[1],ulst[2],ylst[1]-ylst[0],ylst[2]-ylst[0],1)
    tau       = mulst[0] * dudy
    u_tau     = np.sqrt(tau/rholst[0])
    lv       = mulst[0]/rholst[0]/u_tau 
    delta     = 5.2  
    Tinfi     = 160.15 
    print("x = ",line_loc[0])
    print("u_tau = ",u_tau)
    print("tau   = ",tau)
    print("lv    = ",lv)  
    print("rho_w = ",rholst[0])   
    print("==============")                          
    # plot curve
    os.chdir(OutPath)
    out_file = "x_" + str(line_loc[0]) + ".dat"
    with open(out_file,"w") as f:
        for i in range(len(ylst)):
            deltay = ylst[i] - ylst[0]
            yplus  = deltay/lv
            yplus2 = ylst[i]/lv 
            uplus  = ulst[i]/u_tau
            ydelta = ylst[i]/delta
            uunor  = uulst[i]*rholst[i]/rholst[0]/(u_tau*u_tau)
            vvnor  = vvlst[i]*rholst[i]/rholst[0]/(u_tau*u_tau)
            wwnor  = wwlst[i]*rholst[i]/rholst[0]/(u_tau*u_tau)
            uvnor  = uvlst[i]*rholst[i]/rholst[0]/(u_tau*u_tau)
            Tplus  = Tlst[i]/Tinfi
            f.write(str('%.7f'%xlst[i])  +' '+\
                    str('%.7f'%ylst[i])  +' '+\
                    str('%.7f'%deltay)   +' '+\
                    str('%.7f'%yplus)    +' '+\
                    str('%.7f'%yplus2)   +' '+\
                    str('%.7f'%ydelta)   +' '+\
                    str('%.7f'%ulst[i])  +' '+\
                    str('%.7f'%uplus)    +' '+\
                    str('%.7f'%uunor)    +' '+\
                    str('%.7f'%vvnor)    +' '+\
                    str('%.7f'%wwnor)    +' '+\
                    str('%.7f'%uvnor)    +' '+\
                    str('%.7f'%Tlst[i])  +' '+\
                    str('%.7f'%Tplus)    +' '+\
                    str('%.7f'%rholst[i])+' '+\
                    str('%.7f'%uulst[i]) +' '+\
                    str('%.7f'%vvlst[i]) +' '+\
                    str('%.7f'%wwlst[i]) +'\n')
            
#%% Get the x-z plane averaged data within blocks intersected with a line
def GetLineAve(line_loc,zonegrp,FoldPath,OutPath):
    # line_loc = [x1, y1,  x2, y2]
    # zonegrp is a list of instances of Class ZoneGroup
    # set an empty list for storing the intersected zone
    # only vertical line works so far
    xlst   = [line_loc[0]]
    ylst   = [line_loc[1]]
    ulst   = [0.0]
    uulst  = [0.0]
    vvlst  = [0.0]
    wwlst  = [0.0]
    uvlst  = [0.0]
    rholst = [0.0]
    mulst  = [0.0]
    Tlst   = [0.0]
    zone_intersect_i = []
    for i in range(len(zonegrp)):
        line_intersected =  (line_loc[0] >= zonegrp[i].xmin) and \
                            (line_loc[0] <  zonegrp[i].xmax) and \
                            (line_loc[1] <  zonegrp[i].ymax) and \
                            (line_loc[3] >  zonegrp[i].ymin)                        
        if line_intersected:
            # store the index of intersected zones
            zone_intersect_i.append(i)  
    # have selected zones intersected | start to read in specific zones   
    FileList = sorted(get_filelist(FoldPath))
    dataset  = tp.data.load_tecplot_szl(FileList)
    zone     = dataset.zone
    for i in range(len(zone_intersect_i)):
        # For each zone group, x,y DO NOT need to be averaged
        # other variables need to be averaged
        # Firstly get the total number of zones in a zone group
        zone_span_n = len(zonegrp[zone_intersect_i[i]].zonelist)
        for j in range(zone_span_n):
            zonename = zonegrp[zone_intersect_i[i]].zonelist[j] 
            # read in x and y for the first time  
            if j==0:
                Nxi, Nyi, Nzi = dataset.zone(zonename).dimensions
                x      = zone(zonename).values('x').as_numpy_array() 
                y      = zone(zonename).values('y').as_numpy_array() 
                u      = zone(zonename).values('<u>').as_numpy_array()
                uu     = zone(zonename).values('<u`u`>').as_numpy_array()
                vv     = zone(zonename).values('<v`v`>').as_numpy_array()
                ww     = zone(zonename).values('<w`w`>').as_numpy_array()
                uv     = zone(zonename).values('<u`v`>').as_numpy_array()                
                rho    = zone(zonename).values('<rho>').as_numpy_array()
                mu     = zone(zonename).values('<mu>').as_numpy_array()
                T      = zone(zonename).values('<T>').as_numpy_array()
                ugrp   = u
                uugrp  = uu
                vvgrp  = vv
                wwgrp  = ww
                uvgrp  = uv
                rhogrp = rho                
                mugrp  = mu
                Tgrp   = T
            else:
                u      = zone(zonename).values('<u>').as_numpy_array()
                uu     = zone(zonename).values('<u`u`>').as_numpy_array()
                vv     = zone(zonename).values('<v`v`>').as_numpy_array()
                ww     = zone(zonename).values('<w`w`>').as_numpy_array()
                uv     = zone(zonename).values('<u`v`>').as_numpy_array()                
                rho    = zone(zonename).values('<rho>').as_numpy_array()
                mu     = zone(zonename).values('<mu>').as_numpy_array()
                T      = zone(zonename).values('<T>').as_numpy_array()                
                # ugrp stores u for all spanwise zones(u of the group)
                ugrp   = np.vstack((ugrp,u))
                uugrp  = np.vstack((uugrp,uu))
                vvgrp  = np.vstack((vvgrp,vv))
                wwgrp  = np.vstack((wwgrp,ww))
                uvgrp  = np.vstack((uvgrp,uv))             
                rhogrp = np.vstack((rhogrp,rho))
                mugrp  = np.vstack((mugrp,mu))
                Tgrp   = np.vstack((Tgrp,T))
        # spanwise average the variables for zone group i
        # when ugrp has only one dimension, no need to to average
        if np.ndim(ugrp) == 1:
            # if here use np.mean, will collapse to a single value
            umean   = ugrp
            uumean  = uugrp
            vvmean  = vvgrp 
            wwmean  = wwgrp
            uvmean  = uvgrp            
            rhomean = rhogrp
            mumean  = mugrp
            Tmean   = Tgrp         
        else:                     
            umean   = np.mean(ugrp,  axis = 0)
            uumean  = np.mean(uugrp, axis = 0)
            vvmean  = np.mean(vvgrp, axis = 0)
            wwmean  = np.mean(wwgrp, axis = 0)
            uvmean  = np.mean(uvgrp, axis = 0)
            rhomean = np.mean(rhogrp,axis = 0)
            mumean  = np.mean(mugrp, axis = 0)
            Tmean   = np.mean(Tgrp,  axis = 0)
        # reshape the list of variables into matrix    
        x = np.unique(x)
        y = np.unique(y)
        # notice the order of u
        umean   = np.reshape(umean,  (Nyi,Nxi))
        uumean  = np.reshape(uumean, (Nyi,Nxi))
        vvmean  = np.reshape(vvmean, (Nyi,Nxi))
        wwmean  = np.reshape(wwmean, (Nyi,Nxi))
        uvmean  = np.reshape(uvmean, (Nyi,Nxi))
        rhomean = np.reshape(rhomean,(Nyi,Nxi))
        mumean  = np.reshape(mumean, (Nyi,Nxi))
        Tmean   = np.reshape(Tmean,  (Nyi,Nxi))
        # do average along x direction(plane averaged data)    
        if i < (len(zone_intersect_i)-1):
            y_i_max = np.size(y)-1
        else:
            y_i_max = np.size(y)             
        for y_i in range(y_i_max):
            # axis=1 do average in range Nxi, axis=0 in range Nyi
            local_u   = np.mean(umean,  axis=1)
            local_uu  = np.mean(uumean, axis=1)
            local_vv  = np.mean(vvmean, axis=1)
            local_ww  = np.mean(wwmean, axis=1)
            local_uv  = np.mean(uvmean, axis=1)
            local_rho = np.mean(rhomean,axis=1)
            local_mu  = np.mean(mumean, axis=1)            
            local_T   = np.mean(Tmean,  axis=1)
            
            