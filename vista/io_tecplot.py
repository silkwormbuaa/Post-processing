# -*- coding: utf-8 -*-
'''
@File    :   io_tecplot.py
@Time    :   2022/10/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Post processing functions with pytecplot module
'''

import os 
import tecplot           as     tp
import numpy             as     np 
import pandas            as     pd
from   .timer            import timer 
from   .math_opr         import lin_interpo
from   .tools            import get_filelist
from   .tools            import if_overlap_2d
from   .tools            import if_penetrate


# ----------------------------------------------------------------------
# >>> CLASS - ZONE GROUP                                          ( 0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/13  - created
#
# Desc
#
# - A zone group is defined as zones(blocks) sharing the same
#   x-y location
#
# ----------------------------------------------------------------------

class ZoneGroup:

    def __init__( self, xctr, yctr ):
        
        # zone center location
        
        self.xctr = xctr
        self.yctr = yctr
        
        # left bottom and right top corner locations
    
        self.xmin     = 0.0
        self.ymin     = 0.0
        self.xmax     = 0.0
        self.ymax     = 0.0

        # zones' name, which shares the same x-y location
        
        self.zonelist = list()
    
    
# ----------------------------------------------------------------------
# >>> SORTED ZONELIST                                           ( 0.1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/15  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def sort_zonelist( self ):
        
        # check if the zonelist dictionary is empty
        
        if ( len(self.zonelist) == 0 ):
            
            raise ValueError("The zonelist is empty")
        
        else:
            
            self.zonelist = sorted(self.zonelist)


# ----------------------------------------------------------------------
# >>>  GET ZONE GROUP                                             ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/13  - created
#
# Input
# 
# - TP_stat folder
#
# ----------------------------------------------------------------------

def get_zonegrp( folderpath ):
    
    with timer("Read in statistic data"):
        
        # get the all szplt file's path+name   
                 
        filelist = sorted( get_filelist(folderpath) )
        
        print("Number of files in this folder: %s" % len( filelist ))
        
        # read in szplt file using pytecplot
        
        dataset  = tp.data.load_tecplot_szl(filelist)
        zone     = dataset.zone
        
        # zonegrp is a list of instances of class ZoneGroup 
               
        zonegrp  = list()
        
        # zonectr is a list containing the locations of zone centers    
                  
        zonectr  = []
        
        # classify zones one by one
        
        for i in range(np.size(filelist)):
            
            zonename = zone(i).name
            
            # find the range of a zone
            
            xmin = round(dataset.zone(zonename).values('x').min(),7)
            ymin = round(dataset.zone(zonename).values('y').min(),7)
            zmin = round(dataset.zone(zonename).values('z').min(),7)
            xmax = round(dataset.zone(zonename).values('x').max(),7)
            ymax = round(dataset.zone(zonename).values('y').max(),7)
            zmax = round(dataset.zone(zonename).values('z').max(),7)
            xctr = 0.5 * (xmin + xmax)
            yctr = 0.5 * (ymin + ymax)
            zctr = 0.5 * (zmin + zmax)
            
            # if the zone's center is in the list, add this zone's 
            # name to the already existed zone group.
            if [ xctr, yctr ] in zonectr:
                
                indx = zonectr.index([xctr,yctr])
                zonegrp[indx].zonelist.append((zctr,zonename))
                
            # if the zone is not in the list, add a new zone group 
            # to the zone group list, then add this zone to this
            # new zone group    
            
            else:
                
                zonegrp.append(ZoneGroup(xctr,yctr))            
                zonectr.append([xctr,yctr])
                
                indx = zonectr.index([xctr,yctr])
                zonegrp[indx].zonelist.append((zctr,zonename))
                zonegrp[indx].xmin = xmin
                zonegrp[indx].ymin = ymin
                zonegrp[indx].xmax = xmax
                zonegrp[indx].ymax = ymax
                
        # sort the zonelist in each zone group
        
        for group in zonegrp:
            
            group.sort_zonelist()

        # sort the zone center coordinates firstly based on y, then x  
              
        zonectr.sort(key=lambda x: x[1]) 
        
        # sort the zonegrp by location ang output the groups
        
        zonegrp.sort(key=lambda x: (x.yctr,x.xctr))     
                
        # print the zone groups info
                
        os.chdir(folderpath)
        os.chdir(os.pardir)
        
        with open("zonelist_sorted.dat","w") as f:
            
            for group in zonegrp:
                
                print(group.zonelist)
                
                f.write(str('{:<17.8e}'.format(group.xctr)) + ' ')
                f.write(str('{:<17.8e}'.format(group.yctr)) + ' ')
                f.write(str('{:<17.8e}'.format(group.xmin)) + ' ')
                f.write(str('{:<17.8e}'.format(group.ymin)) + ' ')
                f.write(str('{:<17.8e}'.format(group.xmax)) + ' ')
                f.write(str('{:<17.8e}'.format(group.ymax)) + ' ')
                
                for pair in group.zonelist:
                    
                    f.write(str(pair[0]) + ' ')
                
                for pair in group.zonelist:
                    
                    f.write(str(pair[1]) + ' ')
                    
                f.write('\n')    
                
    return zonegrp 


# ----------------------------------------------------------------------
# >>> READ ZONE GROUP LIST                                        ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/13  - created
#
# Input
#
# - zonelist.dat file
# - select mode: all zones or zones with IB
#
# ----------------------------------------------------------------------

def read_zonelist( listfile, wall=False ):
    
    zonegrp = list()
    
    with open( listfile ) as f:
        
        i = 0
        
        for line in f.readlines():
            
            cleanline = line.strip().split()
            
            xctr = float( cleanline[0] )
            yctr = float( cleanline[1] )
            
            if wall is True:
                
                if yctr < 0.0 :
                    
                    zonegrp.append( ZoneGroup(xctr,yctr)  )
                    zonegrp[i].xmin = float( cleanline[2] )
                    zonegrp[i].ymin = float( cleanline[3] )
                    zonegrp[i].xmax = float( cleanline[4] )
                    zonegrp[i].ymax = float( cleanline[5] )
                    
                    for j in range( 6, 14 ):
                        zonegrp[i].zonelist.append( (cleanline[j],cleanline[j+8] ))
                    i = i + 1
                    
            else:
                
                zonegrp.append( ZoneGroup(xctr,yctr)  )
                zonegrp[i].xmin = float( cleanline[2] )
                zonegrp[i].ymin = float( cleanline[3] )
                zonegrp[i].xmax = float( cleanline[4] )
                zonegrp[i].ymax = float( cleanline[5] )
                
                for j in range( 6, 14 ):
                    zonegrp[i].zonelist.append( (cleanline[j],cleanline[j+8] ))
                i = i + 1
                
    return zonegrp


# ----------------------------------------------------------------------
# >>> AVERAGE BLOCK                                              ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/13  - created
#
# Desc
#
# - Average data in a rectangular region, collapse x,z dimension
# - The region right edge does not need to fit block edge
# ----------------------------------------------------------------------

def ave_block( zonegrp, FoldPath, outfile, block_dim=2 ):
    
    RegionRange = [-71.75, 0.0, -60.6875, 22.8044] #22.8044
    
    # Get the zones overlapped with the select region
    
    zone_overlap_indx = []
    layers_yctr = []
    for i in range(len(zonegrp)):
        
        zonegrp_range = [zonegrp[i].xmin, zonegrp[i].ymin,
                         zonegrp[i].xmax, zonegrp[i].ymax]
        
        Overlap = if_overlap_2d(zonegrp_range, RegionRange)

        if Overlap: 
            zone_overlap_indx.append(i)
            layers_yctr.append(zonegrp[i].yctr)
            print(i,[pair[1] for pair in zonegrp[i].zonelist])
            
    # Group the zones by their y location
    
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
         
    # Read in the szplt file with tecplot

    FileList = sorted(get_filelist(FoldPath))
    dataset  = tp.data.load_tecplot_szl(FileList)
    zone     = dataset.zone
    
    # Manipulate the block data, read in variables, reduce dimension
    
    u_ls     = []
    uu_ls    = []
    vv_ls    = []
    ww_ls    = []
    uv_ls    = []
    rho_ls   = []
    T_ls     = []
    y_ls     = []
    
    
    # outer loop to get mean variables along y by combining mean variables 
    # for blocks at different y location.
    
    for i in range(len(zonels_y)):
        
        # loop to get mean variables for blocks at the same y location(x averaging)
        
        for j in range(len(zonels_y[i])):
            
            # inner loop to average variables for blocks at the same x-y location 
                       
            zone_span_n = len(zonegrp[zonels_y[i][j]].zonelist)
            for k in range(zone_span_n):
                zonename = zonegrp[zonels_y[i][j]].zonelist[k][1]
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
                    
            # Span-wise average variables for zone_group i
            
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
                
            # Reshape the list of variables into matrix
            # if block is 3D,need to average block in Z direction firstly
            
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
        
        # for x-stacked zones, find where to slice
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
        # get the x-averaged variables
        umean   = np.mean(uxstack,  axis = 1)
        uumean  = np.mean(uuxstack, axis = 1)
        vvmean  = np.mean(vvxstack, axis = 1)
        wwmean  = np.mean(wwxstack, axis = 1)
        uvmean  = np.mean(uvxstack, axis = 1)
        rhomean = np.mean(rhoxstack,axis = 1)
        Tmean   = np.mean(Txstack,  axis = 1)
#        print(umean)
#        print(umean.shape)
        # Combine the segments along y
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
            # delete the repeted data point
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
    with open(outfile,'w') as f:
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


# ----------------------------------------------------------------------
# >>> TIME AVERAGE f_bl DATA                                      ( 4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/13  - created
#
# Desc
#
# - Do time averaging for f_bl files and
# - return mean values to "statistic_average.dat" file
#
# Input
#
# - folder containing f_bl files, outfile is in pardir
# ----------------------------------------------------------------------

def timeave_fbl( FoldPath ):
    
    # mean(time averaged) variables
    
    fx_m   = []
    px_m   = []
    fx_d_m = []
    px_d_m = []
    rho_m  = []
    mu_m   = []
    nu_m   = []
    surf_m = []
    surf_p = 9.5875
    
    os.chdir( FoldPath )
    FileList = sorted( get_filelist(FoldPath) )
    
    for filename in FileList:
        
        t_b    = []
        fx_b   = []
        fy_b   = []
        fz_b   = []
        px_b   = []
        py_b   = []
        pz_b   = []
        fx_d_b = []
        fy_d_b = []
        fz_d_b = []
        px_d_b = []
        py_d_b = []
        pz_d_b = []
        rho_b  = []
        mu_b   = []
        nu_b   = []
        surf_b = []
        
        with open(filename) as f:
            
            for line in f.readlines():
                
                cleanl = line.strip().split()
                t_b.   append(float(cleanl[0]))
                fx_b.  append(float(cleanl[1]))  
                fy_b.  append(float(cleanl[2]))  
                fz_b.  append(float(cleanl[3]))  
                px_b.  append(float(cleanl[4]))  
                py_b.  append(float(cleanl[5]))  
                pz_b.  append(float(cleanl[6]))  
                fx_d_b.append(float(cleanl[7]))
                fy_d_b.append(float(cleanl[8]))
                fz_d_b.append(float(cleanl[9]))
                px_d_b.append(float(cleanl[10]))
                py_d_b.append(float(cleanl[11]))
                pz_d_b.append(float(cleanl[12]))
                rho_b. append(float(cleanl[13])) 
                mu_b.  append(float(cleanl[14]))  
                nu_b.  append(float(cleanl[15]))  
                surf_b.append(float(cleanl[16]))
                
        fx_m  .append(np.mean(fx_b))
        px_m  .append(np.mean(px_b))
        fx_d_m.append(np.mean(fx_d_b))
        px_d_m.append(np.mean(px_d_b))
        rho_m .append(np.mean(rho_b))
        mu_m  .append(np.mean(mu_b))
        nu_m  .append(np.mean(nu_b))
        surf_m.append(np.mean(surf_b))
        
    # from mean results of each block to get average results        
    
    tau_av1 = (sum(fx_m) + sum(px_m)) / (surf_p*len(surf_m))
    tau_fx1 = sum(fx_m) / (surf_p*len(surf_m))
    tau_px1 = sum(px_m) / (surf_p*len(surf_m))
    rho_av1 = sum(np.multiply(rho_m,surf_m)) / sum(surf_m)
    u_tau1  = np.sqrt(abs(tau_av1) / rho_av1)
    nu_av1  = sum(np.multiply(nu_m,surf_m)) / sum(surf_m)
    mu_av1  = sum(np.multiply(mu_m,surf_m)) / sum(surf_m)
    lv_av1  = nu_av1 / u_tau1
    nu_av2  = mu_av1 / rho_av1
    u_tau2  = np.sqrt((abs(sum(fx_d_m) + sum(px_d_m))) / \
                  (surf_p*len(surf_m)))
    
    # output statistics averaged results to txt file
    
    os.chdir(os.pardir)
    with open("statistic_average.dat",'w') as f:
        f.write('{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}'\
               .format('tau_av1', 'rho_av1', 'u_tau1', 'nu_av1',\
                       'mu_av1' , 'lv_av1' , 'nu_av2', 'u_tau2',\
                       'tau_fx1', 'tau_px1') + '\n')
        f.write(str('{:<17.8e}'.format(tau_av1)))
        f.write(str('{:<17.8e}'.format(rho_av1)))
        f.write(str('{:<17.8e}'.format(u_tau1)))
        f.write(str('{:<17.8e}'.format(nu_av1)))
        f.write(str('{:<17.8e}'.format(mu_av1)))
        f.write(str('{:<17.8e}'.format(lv_av1)))
        f.write(str('{:<17.8e}'.format(nu_av2)))
        f.write(str('{:<17.8e}'.format(u_tau2)))
        f.write(str('{:<17.8e}'.format(tau_fx1)))
        f.write(str('{:<17.8e}'.format(tau_px1)))
        print("finish write statistic results.")


# ----------------------------------------------------------------------
# >>> AVERAGE BLOCK WITH IB                                       ( 5 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/15  - created
#
# Desc
#
# - average blocks into a vertical line for plotting profile 
# - like ave_block(), but can drop values beneath wall
#
# Input
#
# - spf_ave : if do superficial averaging, otherwise arithmatic averaging
# ----------------------------------------------------------------------

def ave_block_ib( folderpath, zonegrps, 
                  outfile,    regionrange,
                  var_list=None,
                  period_ave=None,
                  zone_row=None,
                  spf_ave=False):
    
    # get the zones overlapped with the select region
    
    zone_overlap_indx = []
#    layers_yctr = []
    
    for i in range( len(zonegrps) ):
        
        zonegrp_range = [zonegrps[i].xmin, zonegrps[i].ymin,
                         zonegrps[i].xmax, zonegrps[i].ymax]
        
        Overlap = if_overlap_2d( zonegrp_range, regionrange )

        if Overlap: 
            zone_overlap_indx.append(i)
#            layers_yctr.append(zonegrp[i].yctr)
            print(i,[pair[1] for pair in zonegrps[i].zonelist])
    
      
#    layer_n = len( np.unique(layers_yctr))
    zonels_y = []
    
    for i in range(7):
        zonels_y.append([])
        
    for i in range(len(zone_overlap_indx)):
        if (zonegrps[zone_overlap_indx[i]].yctr < 0.0):
            zonels_y[0].append(zone_overlap_indx[i])
        elif (0.0        < zonegrps[zone_overlap_indx[i]].yctr < 1.73695342):
            zonels_y[1].append(zone_overlap_indx[i])
        elif (1.73695342 < zonegrps[zone_overlap_indx[i]].yctr < 5.01031265):
            zonels_y[2].append(zone_overlap_indx[i])        
        elif (5.01031265 < zonegrps[zone_overlap_indx[i]].yctr < 11.1790910):
            zonels_y[3].append(zone_overlap_indx[i])
        elif (11.1790910 < zonegrps[zone_overlap_indx[i]].yctr < 22.8044042):
            zonels_y[4].append(zone_overlap_indx[i])   
        elif (22.8044042 < zonegrps[zone_overlap_indx[i]].yctr < 44.7127788):
            zonels_y[5].append(zone_overlap_indx[i])
        elif (44.7127788 < zonegrps[zone_overlap_indx[i]].yctr < 86.0000000):
            zonels_y[6].append(zone_overlap_indx[i])
   
    print(zonels_y)
       
    
    filelist = sorted( get_filelist(folderpath) )
    dataset  = tp.data.load_tecplot_szl( filelist )
                
    if var_list is None:
        var_list = [v.name for v in dataset.variables()]
    

    
    for i in range( len(zonels_y) ):

        for j in range(len(zonels_y[i])):
            
            zone_span_n = len(zonegrps[zonels_y[i][j]].zonelist)
            
#            check_zone = []
            
            for k in range(zone_span_n):
                zonename = zonegrps[zonels_y[i][j]].zonelist[k][1]
                
#                check_zone.append(zonename)
                
                for m in range(np.size(var_list)):
                    
                    var = dataset.variable(var_list[m])
                    
                    if m == 0:
                        var_col = var.values(zonename).as_numpy_array()
                    else:
                        var_index = var.values(zonename).as_numpy_array()
                        var_col = np.column_stack((var_col,var_index))
                if zone_row is None:
#                    if ( np.size(dataset.solution_times) ==0 ):
#                        sol_time = 0.0
#                    else:
#                        sol_time = dataset.solution_times[0]
                    zone_row = var_col
                else:
                    zone_row = np.row_stack((zone_row, var_col))
#                    print(len(var_col))
#            print(check_zone)
#            print(len(zone_row))
    df = pd.DataFrame( data=zone_row, columns=var_list )   

    print(len(df.index))
    
    # drop overlapped points and only keep the first 
    
    df = df.groupby(by=['x','y','z'],as_index=False).first()
    
    print(df)
    
    # set values for walldist < 0 as 0.0 for superficial averaging
    
    if spf_ave is True:
        df['<u>'].loc[df['walldist']<0] = 0.0
        df['<u`u`>'].loc[df['walldist']<0] = 0.0
        df['<v`v`>'].loc[df['walldist']<0] = 0.0
        df['<w`w`>'].loc[df['walldist']<0] = 0.0
        df['<u`v`>'].loc[df['walldist']<0] = 0.0
        df['<rho>'].loc[df['walldist']<0] = 0.0
        df['<T>'].loc[df['walldist']<0] = 0.0
    
    # else drop values for walldist < 0 for arithmatic averaging
    
    else:   
        df = df.drop(df[ df['walldist']<0.0 ].index)
    
    print(df)
    
    # drop values that are out of the range
    
    df = df.drop(df[ (df['x'] < regionrange[0]) | 
                     (df['x'] > regionrange[2])   ].index)
    
    print(df)
    
    # averaging over x-z plane
    
    df = df.groupby( by=['y'],as_index=False).mean()
    
    print(df)
    
    # write results to mean_result.dat
    
    os.chdir(folderpath)
    os.chdir(os.pardir)
    
    with open(outfile, 'w') as f:
        
        for var_name in df.columns:
            f.write('{:<17}'.format(var_name))
            
        f.write('\n')
        
        for i in range(df.shape[0]):
            for var_name in df.columns:
                f.write(str('{:<17.8e}'.format(df.at[i,var_name])))
            f.write('\n')
    
    return df       


# ----------------------------------------------------------------------
# >>> SPANWISE PERIODIC AVERAGING                                ( 6 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/18  - created
#
# Desc
#
# ----------------------------------------------------------------------

def spanwise_periodic_ave( dataset,
                           zonegrp,
                           period,
                           width_out,
                           var_list=None,
                           zone_row=None  ):
    
    # loop over all zones (8) within one zonegrp
    for zone in zonegrp.zonelist:
        
        # check if output variables are specified
        
        if var_list is None:
            
            var_list = [ v.name for v in dataset.variables() ]
        
        # read in the variables with tecplot dataset
        # first read in and stack variable columns within one zone
        
        for i in range( np.size(var_list) ):
            
            var = dataset.variable( var_list[i] )
            
            if i == 0:
                
                # zone[1] is 'B00xxxx', zone[0] is z location of zone center
                
                var_col = var.values( zone[1] ).as_numpy_array()
            
            else:
                
                var_index = var.values( zone[1] ).as_numpy_array()
                var_col = np.column_stack( (var_col,var_index) )

        # then stack zones(vertically)

        if zone_row is None:
            
            zone_row = var_col
            
        else:
            
            zone_row = np.row_stack( (zone_row, var_col) )
    
    # store data into pandas
    
    df = pd.DataFrame( data=zone_row, columns=var_list ) 

    # drop overlapped data at blocks' interface

    df = df.drop_duplicates( subset = ['x','y','z'], keep='first' )
    
    x = np.array( df['x'] )
    y = np.array( df['y'] )
    z = np.array( df['z'] )
    
    Nx = len( np.unique(x) )
    Ny = len( np.unique(y) )
    Nz = len( np.unique(z) )

    # doing spanwise periodic averaging, need to manipulate z
    # first, reshape z into a 3D matrix
    
    z = np.reshape( z, (Nz,Ny,Nx) )
    
    # get the width (length in z) of this zone group
    
    width = max(np.unique(z)) - min(np.unique(z))

    # number of periods over the whole span
    
    n_period = round(width/period)
    
    # how many points are in one period

    n_zpoints = round( (Nz-1)/n_period )
    
    # z_frame is the new z coordinate within one period
    # z_frame need to drop the last points to avoid repetitive 
    # writing of period interface data
    
    z_frame = np.linspace( 0.0, period, n_zpoints+1)
    z_frame = z_frame[ 0:-1 ]
    z_frame = z_frame.round( decimals=5 )   # maybe useless
    
    # assigning x-y planes with new repetitive z coordinates
    
    for i in range( len(z)-1 ):
        
        index = i % n_zpoints
        
        z[i,:,:] = z_frame[ index ]
        
    z[-1,:,:] = z_frame[0]    

    # flatten z back into a 1D array, for covenience of being added to pandas
    
    z = z.ravel()
    
    df['z'] = z
    
    # finish manipulating z coordinates, z range only within one periodic
    # now can do grouping and averaging
    
    df = df.groupby( by=['z','y','x'], as_index=False ).mean()
    
    # number of peroid for output file(depending on how width of out zone)
     
    N_period = round( width_out/period )
    
    # copy dataframe to help assignvalues for periodic zones
    
    df_cp   = df.copy( deep=True )
    z_cp    = np.array( df['z'] )
    
    # leave a copy of the first x-y plane, later will be added to tail
    
    df_tail = df.iloc[ 0 : Nx*Ny ].copy()
    
    # assign z coordiantes for the tail x-y layer
    # (the only difference between the first and last x-y layer)
    
    z_tail = np.array( [width_out] * (Nx*Ny) )
    
    df_tail['z'] = z_tail
    
    # when output file requires not only one period
    # need to write extra repetitive; N always >= 1, no worry
    
    for i in range( 1, N_period ):
        
        # extra block need to shif z 
        
        z_extra = np.add( z_cp, i*period )
        
        df_cp['z'] = z_extra 
        
        df = pd.concat([df,df_cp])
    
    # get final dataframe

    df = pd.concat([df,df_tail])
    
    # adjust index and order of columns variables before output
    
    df = df.reset_index(drop=True)  
    df = df[var_list]
    
    print(df)
    
    return df


# ----------------------------------------------------------------------
# >>> Streamwise Mean Line                                       ( 7 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/24  - created
#
# Desc
#
# - get a streamwise line with spanwise mean results at certain y_0
#
# ----------------------------------------------------------------------

def get_xline( folderpath, zonegrps,
               segment, 
               outfile,
               var_list=None,
               spf_ave=False) : 
    
    # values used for normalization
    
    p_i   = 45447.289
    x_imp = 50.4
    delta = 5.2
    
    # get the zones penetrated by the segment
    
    zone_penetr_indx = []
    zonels = []
    
    for i in range(len(zonegrps)):
        
        zonegrp_range = [zonegrps[i].xmin, zonegrps[i].ymin,
                         zonegrps[i].xmax, zonegrps[i].ymax]
        
        Penetrate = if_penetrate( zonegrp_range, segment )
        
        if Penetrate:
            zone_penetr_indx.append(i)
            zonels.append( [pair[1] for pair in zonegrps[i].zonelist] )
            print(i,zonegrps[i].xctr,[pair[1] for pair in zonegrps[i].zonelist])
    
#    print(zonels)

    filelist = sorted( get_filelist(folderpath) )
    dataset  = tp.data.load_tecplot_szl( filelist )
                
    if var_list is None:
        var_list = [v.name for v in dataset.variables()]
        
    zone_row = None
    
    for i in range( (len(zonels)) ):
        
        zone_span_n = len(zonels[i])
        
        for j in range( zone_span_n ):
            
            zonename = zonels[i][j]
            
            for k in range( np.size(var_list) ):
                
                var = dataset.variable(var_list[k])
                
                if k == 0:
                    var_col = var.values( zonename ).as_numpy_array()
                else:
                    var_index = var.values( zonename ).as_numpy_array()
                    var_col = np.column_stack(( var_col, var_index))
            
            if zone_row is None:
                
                zone_row = var_col
            
            else:
                
                zone_row = np.row_stack(( zone_row, var_col ))
    
    df = pd.DataFrame( data=zone_row, columns=var_list )
    
#    print(len(df.index))
    
    # drop repetitive data
    
    df = df.groupby(by=['x','y','z'], as_index=False).first()
    
    # drop points in solid
    
    if spf_ave is True:
        df['p`p`'].loc[df['walldist']<0] = 0.0
    else:
        df = df.drop(df[ df['walldist']<0.0 ].index)
    
    # do spanwise averaging
    
    df = df.groupby( by=['x','y'], as_index=False).mean()
    
    # get the interval where the target y falls
    
    y = np.array( df['y'] )
    
    y = np.sort( np.unique(y) )
    
    print(y)
    
    if (segment[0] < min(y)) or (segment[0] > max(y)):
        
        raise ValueError("y_loc is out of data's range")
    
    else:
        
        for index in range( len(y) ):
            if (y[index] >= segment[0]):
                break
    
    print(index)
    
    df1 = df.loc[df['y']==y[index-1]].copy()
    df2 = df.loc[df['y']==y[index]].copy()
    
    df1['<p`>'] = df1['<p`p`>'].apply(np.sqrt)
    df2['<p`>'] = df2['<p`p`>'].apply(np.sqrt)
    
    p1 = np.array(df1['<p`>'])
    p2 = np.array(df2['<p`>'])
    
    p_interpo = []
    
    for i in range(len(p1)):
        p = lin_interpo(segment[0],y[index-1],p1[i],y[index],p2[i])
        p_interpo.append(p)
    
    x = np.array(df1['x'])
    df_out = pd.DataFrame(data=x,columns=['x'])
    df_out['<p`>'] = p_interpo
    
    # normalization
    
    df_out['(x-x_imp)/δ'] = np.subtract(x,x_imp)/delta
    
    df_out['<p`>_'] = np.divide(p_interpo,p_i)

    os.chdir(folderpath)
    os.chdir(os.pardir)
    
    fmt = '%.7f'
    df_out.to_csv(outfile, sep=' ', float_format=fmt, index=False)
        
    
    print(df_out)
#    y = np.array(df['y'])
#    
#    y = np.unique(y)
#    
#    print(y)
    