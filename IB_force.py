# -*- coding: utf-8 -*-
'''
@File    :   get_IB.py
@Time    :   2022/09/23 15:04:47
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os

import numpy             as     np


#%% Read in forces, pressure... on immersed boundary.
def ReadAForce(filename):
#---mean(time averaged) variables    
    fx_m   = 0.0
    px_m   = 0.0
#    fx_d_m = 0.0
#    px_d_m = 0.0
#    rho_m  = 0.0
#    mu_m   = 0.0
#    nu_m   = 0.0
#    surf_m = 0.0
    #surface of plane
    surf_p = 9.5875
#---containers for read in data
    t_b    = []
    fx_b   = []
    fy_b   = []
    fz_b   = []
    px_b   = []
    py_b   = []
    pz_b   = []
#    fx_d_b = []
#    fy_d_b = []
#    fz_d_b = []
#    px_d_b = []
#    py_d_b = []
#    pz_d_b = []
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
#                fx_d_b.append(float(cleanl[7]))
#                fy_d_b.append(float(cleanl[8]))
#                fz_d_b.append(float(cleanl[9]))
#                px_d_b.append(float(cleanl[10]))
#                py_d_b.append(float(cleanl[11]))
#                pz_d_b.append(float(cleanl[12]))
                rho_b. append(float(cleanl[13])) 
                mu_b.  append(float(cleanl[14]))  
                nu_b.  append(float(cleanl[15]))  
                surf_b.append(float(cleanl[16]))
    fx_m   = np.mean(fx_b)
    px_m   = np.mean(px_b)
    py_m   = np.mean(py_b)
#    fx_d_m = np.mean(fx_d_b)
#    px_d_m = np.mean(px_d_b)
    rho_m  = np.mean(rho_b)
    mu_m   = np.mean(mu_b)
    nu_m   = np.mean(nu_b)
    surf_m = np.mean(surf_b)
#---get tau and 
    tau_w  = (fx_m + px_m)/surf_p
    p_w    = py_m/surf_p
    
    return tau_w, p_w

#%% def a method to get points value of Cf,Pf
def GetIBForce(zonegrp,ForceFolderPath,Outfile):

    d_p = 0.5*0.9886*507*507
    p_i = 45447.289
    x   = []
    fx  = []
    p   = []
    for zones in zonegrp:
        x.append(zones.xctr)
#        print(zones.zonelist)
        fx_local = []
        p_local  = []
        for blockname in zones.zonelist:
            filename = "f_bl_" + blockname.lstrip("B") + ".dat"
            tau_w, p_w = ReadAForce(ForceFolderPath+filename)
            fx_local.append(tau_w)
            p_local.append(p_w)
        fx.append(np.mean(fx_local))
        p .append(np.mean(p_local))
    Cf = -np.divide(fx,d_p)
    Pf = np.divide(p,p_i)
    
    with open(Outfile,'w') as f:
        f.write("x           Cf        Cp\n")
        for i in range(len(x)):
            f.write(str('%.7f'%x[i])+' '+str('%.7f'%Cf[i])+' '\
                + str('%.7f'%Pf[i])+'\n')

