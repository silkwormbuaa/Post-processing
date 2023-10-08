# -*- coding: utf-8 -*-
'''
@File    :   ReadIn.py
@Time    :   2022/07/06 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Class and functions related to read in data
'''

import os

from   os                import path

import tecplot           as     tp

import numpy             as     np

from   vista.timer       import timer

from   vista.tools       import get_filelist
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

#%%  get the zone groups(same x-y location) from statistic data.
def GetZonegrp(FoldPath):
    with timer("Read in statistic data"):
        # get the all szplt file's path+name            
        FileList = sorted(get_filelist(FoldPath))
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
def ReadZonegrpWall(filename):
    zonegrp = list()
    with open(filename) as f:
        i = 0
        for line in f.readlines():
            cleanline = line.strip().split()
            xctr = float(cleanline[0])
            yctr = float(cleanline[1])
            if yctr < 0.0 :
                zonegrp.append(ZoneGroup(xctr,yctr))
                zonegrp[i].xmin = float(cleanline[2])
                zonegrp[i].ymin = float(cleanline[3])
                zonegrp[i].xmax = float(cleanline[4])
                zonegrp[i].ymax = float(cleanline[5])
                for j in range(6,14):
                    zonegrp[i].zonelist.append(cleanline[j])
                i = i + 1
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
def ReadBlock(zonegrp,FoldPath,filename,block_dim=2):
    RegionRange = [-71.75, 0.0, -60.6875, 22.8044] #22.8044
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
    FileList = sorted(get_filelist(FoldPath))
    dataset  = tp.data.load_tecplot_szl(FileList)
    zone     = dataset.zone
    
#---Manipulate the block data, read in variables, reduce dimension
    u_ls     = []
    uu_ls    = []
    vv_ls    = []
    ww_ls    = []
    uv_ls    = []
    rho_ls   = []
    T_ls     = []
    y_ls     = []
#---outer loop to get mean variables along y by combining mean variables 
#   for blocks at different y location.
    for i in range(len(zonels_y)):
#---loop to get mean variables for blocks at the same y location(x averaging)
        for j in range(len(zonels_y[i])):
#---inner loop to average variables for blocks at the same x-y location            
            zone_span_n = len(zonegrp[zonels_y[i][j]].zonelist)
            for k in range(zone_span_n):
                zonename = zonegrp[zonels_y[i][j]].zonelist[k]
                u      = zone(zonename).values('<u>').as_numpy_array()
                uu     = zone(zonename).values('<u`u`>').as_numpy_array()
                vv     = zone(zonename).values('<v`v`>').as_numpy_array()
                ww     = zone(zonename).values('<w`w`>').as_numpy_array()
                uv     = zone(zonename).values('<u`v`>').as_numpy_array()
                rho    = zone(zonename).values('<rho>').as_numpy_array()
                T      = zone(zonename).values('<T>').as_numpy_array()
                if k==0:
                    Nxi, Nyi, Nzi = dataset.zone(zonename).dimensions
                    x      = zone(zonename).values('x').as_numpy_array() 
                    y      = zone(zonename).values('y').as_numpy_array() 
                    ugrp   = u
                    uugrp  = uu
                    vvgrp  = vv
                    wwgrp  = ww
                    uvgrp  = uv
                    rhogrp = rho
                    Tgrp   = T
                else:
                    ugrp   = np.vstack((ugrp,  u))
                    uugrp  = np.vstack((uugrp, uu))
                    vvgrp  = np.vstack((vvgrp, vv))
                    wwgrp  = np.vstack((wwgrp, ww))
                    uvgrp  = np.vstack((uvgrp, uv))
                    rhogrp = np.vstack((rhogrp,rho))
                    Tgrp   = np.vstack((Tgrp,  T))
#---Span-wise average variables for zone_group i
            if np.ndim(ugrp) == 1:
                umean_z      = ugrp
                uumean_z     = uugrp
                vvmean_z     = vvgrp
                wwmean_z     = wwgrp
                uvmean_z     = uvgrp
                rhomean_z    = rhogrp
                Tmean_z      = Tgrp
            else:
                umean_z      = np.mean(ugrp,  axis = 0)
                uumean_z     = np.mean(uugrp, axis = 0)
                vvmean_z     = np.mean(vvgrp, axis = 0)
                wwmean_z     = np.mean(wwgrp, axis = 0)
                uvmean_z     = np.mean(uvgrp, axis = 0)
                rhomean_z    = np.mean(rhogrp,axis = 0)
                Tmean_z      = np.mean(Tgrp, axis = 0)
#---Reshape the list of variables into matrix
#--- if block is 3D,need to average block in Z direction firstly
            if block_dim ==3 :
                x = np.reshape(x,(Nzi,Nyi*Nxi))
                y = np.reshape(y,(Nzi,Nyi*Nxi))
                umean_z   = np.reshape(umean_z,  (Nzi,Nyi*Nxi))
                uumean_z  = np.reshape(uumean_z, (Nzi,Nyi*Nxi))
                vvmean_z  = np.reshape(vvmean_z, (Nzi,Nyi*Nxi))
                wwmean_z  = np.reshape(wwmean_z, (Nzi,Nyi*Nxi))
                uvmean_z  = np.reshape(uvmean_z, (Nzi,Nyi*Nxi))
                rhomean_z = np.reshape(rhomean_z,(Nzi,Nyi*Nxi))
                Tmean_z   = np.reshape(Tmean_z,  (Nzi,Nyi*Nxi))
                umean_z   = np.mean(umean_z,  axis=0)
                uumean_z  = np.mean(uumean_z, axis=0)
                vvmean_z  = np.mean(vvmean_z, axis=0)
                wwmean_z  = np.mean(wwmean_z, axis=0)
                uvmean_z  = np.mean(uvmean_z, axis=0)
                rhomean_z = np.mean(rhomean_z,axis=0)
                Tmean_z   = np.mean(Tmean_z,  axis=0)
            x = np.unique(x)
            y = np.unique(y)
            umean_z   = np.reshape(umean_z,  (Nyi,Nxi))
            uumean_z  = np.reshape(uumean_z, (Nyi,Nxi))
            vvmean_z  = np.reshape(vvmean_z, (Nyi,Nxi))
            wwmean_z  = np.reshape(wwmean_z, (Nyi,Nxi))
            uvmean_z  = np.reshape(uvmean_z, (Nyi,Nxi))
            rhomean_z = np.reshape(rhomean_z,(Nyi,Nxi))
            Tmean_z   = np.reshape(Tmean_z, (Nyi,Nxi))
            if j == 0:
                uxstack   = umean_z
                uuxstack  = uumean_z
                vvxstack  = vvmean_z
                wwxstack  = wwmean_z
                uvxstack  = uvmean_z
                rhoxstack = rhomean_z
                Txstack   = Tmean_z
                xxstack   = x
            else:
                #delete the last element of previous block
                uxstack  =np.hstack((np.delete(uxstack,  -1,axis=1),umean_z))
                uuxstack =np.hstack((np.delete(uuxstack, -1,axis=1),uumean_z))
                vvxstack =np.hstack((np.delete(vvxstack, -1,axis=1),vvmean_z))
                wwxstack =np.hstack((np.delete(wwxstack, -1,axis=1),wwmean_z))
                uvxstack =np.hstack((np.delete(uvxstack, -1,axis=1),uvmean_z))
                rhoxstack=np.hstack((np.delete(rhoxstack,-1,axis=1),rhomean_z))
                Txstack = np.hstack((np.delete(Txstack,  -1,axis=1),Tmean_z))
                xxstack = np.hstack((np.delete(xxstack,  -1),x))
#---for x-stacked zones, find where to slice
        for ii in range(len(xxstack)):
            if xxstack[ii] > RegionRange[2]:
                break         
            if ii == (len(xxstack) - 1):
                ii = len(xxstack)
#        print(uxstack.shape)
#        print(xxstack.shape)
        uxstack   = uxstack[:,:ii]
        uuxstack  = uuxstack[:,:ii]
        vvxstack  = vvxstack[:,:ii]
        wwxstack  = wwxstack[:,:ii]
        uvxstack  = uvxstack[:,:ii]
        rhoxstack = rhoxstack[:,:ii]
        Txstack   = Txstack[:,:ii]
        xxstack   = xxstack[:ii]
#        print(xxstack)
#        print(xxstack.shape)
#---get the x-averaged variables
        umean   = np.mean(uxstack,  axis = 1)
        uumean  = np.mean(uuxstack, axis = 1)
        vvmean  = np.mean(vvxstack, axis = 1)
        wwmean  = np.mean(wwxstack, axis = 1)
        uvmean  = np.mean(uvxstack, axis = 1)
        rhomean = np.mean(rhoxstack,axis = 1)
        Tmean   = np.mean(Txstack,  axis = 1)
#        print(umean)
#        print(umean.shape)
#---Combine the variables along y
        if len(u_ls) == 0:
            u_ls   = np.concatenate((u_ls,  umean))
            uu_ls  = np.concatenate((uu_ls, uumean))
            vv_ls  = np.concatenate((vv_ls, vvmean))
            ww_ls  = np.concatenate((ww_ls, wwmean))
            uv_ls  = np.concatenate((uv_ls, uvmean))
            rho_ls = np.concatenate((rho_ls,rhomean))
            T_ls   = np.concatenate((T_ls,  Tmean))
            y_ls   = np.concatenate((y_ls,  y))
        else:
            u_ls   = np.concatenate((np.delete(u_ls,  -1),  umean))
            uu_ls  = np.concatenate((np.delete(uu_ls, -1), uumean))
            vv_ls  = np.concatenate((np.delete(vv_ls, -1), vvmean))
            ww_ls  = np.concatenate((np.delete(ww_ls, -1), wwmean))
            uv_ls  = np.concatenate((np.delete(uv_ls, -1), uvmean))
            rho_ls = np.concatenate((np.delete(rho_ls,-1), rhomean))
            T_ls   = np.concatenate((np.delete(T_ls,  -1), Tmean))
            y_ls   = np.concatenate((np.delete(y_ls,  -1), y))
    os.chdir(FoldPath)
    os.chdir(os.pardir)
    with open(filename,'w') as f:
        f.write('{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}'\
               .format('y', 'u',  'u`u`','v`v`','w`w`','u`v`',\
               'rho','T') + '\n')
        for i in range(len(y_ls)):
            f.write(str('{:<17.8e}'.format(y_ls[i]))  )
            f.write(str('{:<17.8e}'.format(u_ls[i]))  )
            f.write(str('{:<17.8e}'.format(uu_ls[i])) )
            f.write(str('{:<17.8e}'.format(vv_ls[i])) )
            f.write(str('{:<17.8e}'.format(ww_ls[i])) )
            f.write(str('{:<17.8e}'.format(uv_ls[i])) )
            f.write(str('{:<17.8e}'.format(rho_ls[i])))
            f.write(str('{:<17.8e}'.format(T_ls[i]))  + '\n')

