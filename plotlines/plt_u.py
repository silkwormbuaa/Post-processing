# -*- coding: utf-8 -*-
'''
@File    :   plt_u.py
@Time    :   2022/08/15 22:59:19
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   plot the velocity profile 
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

#%% Read in data first
FoldPath = "/home/wencanwu/my_simulation/temp/220825_lowRe/"
OutPath  = "/home/wencanwu/my_simulation/temp/220825_lowRe/DataPost/"
ForceFoldPath = "/home/wencanwu/my_simulation/temp/220825_lowRe/forces_oneblock/forces_3"
FlatFolder = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/DataPost"
#FoldPath2 = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/"
#ForceFoldPath2 = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/forces/forces_1408"

plt_u   = True
plt_RS  = True
plt_T   = False
Compare = True
#---Read in averaged flow data
os.chdir(FoldPath)
y_ls   = []
u_ls   = []
uu_ls  = []
vv_ls  = []
ww_ls  = []
uv_ls  = []
rho_ls = []
T_ls   = []
with open("mean_result.dat") as f:
    #skip the first line    
    line = f.readline()
    line = f.readline()
    while line:
        cleanl = line.strip().split()
        y_ls.  append(float(cleanl[0]))
        u_ls.  append(float(cleanl[1]))
        uu_ls. append(float(cleanl[2]))
        vv_ls. append(float(cleanl[3]))
        ww_ls. append(float(cleanl[4]))
        uv_ls. append(float(cleanl[5]))
        rho_ls.append(float(cleanl[6]))
        T_ls.  append(float(cleanl[7]))
        #read next line until end        
        line = f.readline()
#---Read in averaged statistic data for normalization
os.chdir(ForceFoldPath)
os.chdir(os.pardir)
with open("statistic_average.dat") as f:
    line    = f.readline()
    line    = f.readline()
    cleanl  = line.strip().split()
    tau_av1 = float(cleanl[0])
    rho_av1 = float(cleanl[1])
    u_tau1  = float(cleanl[2])
    nu_av1  = float(cleanl[3])
    mu_av1  = float(cleanl[4])
    lv_av1  = float(cleanl[5])
    nu_av2  = float(cleanl[6])
    u_tau2  = float(cleanl[7])
#---Get normalized variables
T_inf   = 160.15
delta   = 5.2
yplus   = np.array(y_ls)/lv_av1
uplus   = np.array(u_ls)/u_tau1
uuplus  = np.multiply(uu_ls,rho_ls)/abs(tau_av1)
vvplus  = np.multiply(vv_ls,rho_ls)/abs(tau_av1)
wwplus  = np.multiply(ww_ls,rho_ls)/abs(tau_av1)
uvplus  = np.multiply(uv_ls,rho_ls)/abs(tau_av1)
ydelta  = np.array(y_ls)/delta
Tnorm   = np.array(T_ls)/T_inf
#%%---read spanwise wavy wall case
'''
os.chdir(FoldPath2)
spw_y_ls   = []
spw_u_ls   = []
spw_uu_ls  = []
spw_vv_ls  = []
spw_ww_ls  = []
spw_uv_ls  = []
spw_rho_ls = []
spw_T_ls   = []
with open("mean_result.dat") as f:
    #skip the first line    
    line = f.readline()
    line = f.readline()
    while line:
        cleanl = line.strip().split()
        spw_y_ls.  append(float(cleanl[0]))
        spw_u_ls.  append(float(cleanl[1]))
        spw_uu_ls. append(float(cleanl[2]))
        spw_vv_ls. append(float(cleanl[3]))
        spw_ww_ls. append(float(cleanl[4]))
        spw_uv_ls. append(float(cleanl[5]))
        spw_rho_ls.append(float(cleanl[6]))
        spw_T_ls.  append(float(cleanl[7]))
        #read next line until end        
        line = f.readline()
#---Read in averaged statistic data for normalization
os.chdir(ForceFoldPath2)
os.chdir(os.pardir)
with open("statistic_average.dat") as f:
    line    = f.readline()
    line    = f.readline()
    cleanl  = line.strip().split()
    tau_av1 = float(cleanl[0])
    rho_av1 = float(cleanl[1])
    u_tau1  = float(cleanl[2])
    nu_av1  = float(cleanl[3])
    mu_av1  = float(cleanl[4])
    lv_av1  = float(cleanl[5])
    nu_av2  = float(cleanl[6])
    u_tau2  = float(cleanl[7])
#---Get normalized variables
T_inf   = 160.15
delta   = 5.2
spw_yplus   = np.array(spw_y_ls)/lv_av1
spw_uplus   = np.array(spw_u_ls)/u_tau1
spw_uuplus  = np.multiply(spw_uu_ls,spw_rho_ls)/abs(tau_av1)
spw_vvplus  = np.multiply(spw_vv_ls,spw_rho_ls)/abs(tau_av1)
spw_wwplus  = np.multiply(spw_ww_ls,spw_rho_ls)/abs(tau_av1)
spw_uvplus  = np.multiply(spw_uv_ls,spw_rho_ls)/abs(tau_av1)
spw_ydelta  = np.array(spw_y_ls)/delta
spw_Tnorm   = np.array(spw_T_ls)/T_inf
'''
#%%---Read in flat plate result
if Compare:
    os.chdir(FlatFolder)
    yplus_f   = []
    uplus_f   = []
    uuplus_f  = []
    vvplus_f  = []
    wwplus_f  = []
    uvplus_f  = []
    ydelta_f  = []
    Tnorm_f   = []
    with open("x_-68.0625.dat") as f:
        line = f.readline()
        while line:
            cleanl = line.strip().split()
            yplus_f .append(float(cleanl[4]))
            ydelta_f.append(float(cleanl[5]))
            uplus_f .append(float(cleanl[7]))
            uuplus_f.append(float(cleanl[8]))
            vvplus_f.append(float(cleanl[9]))
            wwplus_f.append(float(cleanl[10]))
            uvplus_f.append(float(cleanl[11]))
            Tnorm_f .append(float(cleanl[13]))
            #read nex line until end
            line = f.readline()

#%% plot u profile
os.chdir(OutPath)
if plt_u :
    fig, ax = plt.subplots(figsize=[10,8])
    if Compare:
        ax.plot(yplus_f,uplus_f,'gray',label=r'$smooth$',ls='--')
        ax.plot(yplus,  uplus, 'blue', label=r'$D=0.5\delta_0$', ls='-')
#        ax.plot(spw_yplus,  spw_uplus, 'green', label=r'$spanwise\ wave$', ls='-')
    else:
        ax.plot(yplus,uplus,label=r'$u^+$')
    
    ax.minorticks_on()
    ax.set_xscale("symlog",linthresh = 1)
    ax.set_xlabel("$y^+$",fontdict={'size':24})  
    ax.tick_params(axis='x',labelsize=15)
    
    ax.set_ylabel(r'$u^+$',fontdict={'size':24})
    ax.tick_params(axis='y',labelsize=15)
    ax.set_xlim([1,1000])
    x_minor = matplotlib.ticker.LogLocator(base=10.0, subs = np.arange(1.0,10.0))
    ax.xaxis.set_minor_locator(x_minor)

    ax.legend(prop={'size':20}) 
    ax.set_title(r"$u^+$ profile",size=20)

    ax.grid()

    plt.savefig("velocity_profile_wavywall_plusspanwise")
    plt.show()

#%% plot Reynolds Stress profile 
if plt_RS :
    fig, ax = plt.subplots(figsize=[10,8])
    if Compare:
        ax.plot(yplus_f,uuplus_f,'b',label=r'$u^\prime u^\prime$ smooth',ls="-")
        ax.plot(yplus_f,uvplus_f,'y',label=r'$u^\prime v^\prime$ smooth',ls="-")
        ax.plot(yplus_f,vvplus_f,'g',label=r'$v^\prime v^\prime$ smooth',ls="-")
        ax.plot(yplus_f,wwplus_f,'r',label=r'$w^\prime w^\prime$ smooth',ls="-")        
        ax.plot(yplus,uuplus,'b',label=r'$u^\prime u^\prime \ D=0.5\delta_0$',ls="--")
        ax.plot(yplus,uvplus,'y',label=r'$u^\prime v^\prime \ D=0.5\delta_0$',ls="--")
        ax.plot(yplus,vvplus,'g',label=r'$v^\prime v^\prime \ D=0.5\delta_0$',ls="--")
        ax.plot(yplus,wwplus,'r',label=r'$w^\prime w^\prime \ D=0.5\delta_0$',ls="--")
#        ax.plot(spw_yplus,spw_uuplus,'b',label=r'$u^\prime u^\prime spanwise\ wave$',ls=":")
#        ax.plot(spw_yplus,spw_uvplus,'y',label=r'$u^\prime v^\prime spanwise\ wave$',ls=":")
#        ax.plot(spw_yplus,spw_vvplus,'g',label=r'$v^\prime v^\prime spanwise\ wave$',ls=":")
#        ax.plot(spw_yplus,spw_wwplus,'r',label=r'$w^\prime w^\prime spanwise\ wave$',ls=":")
    else:    
        ax.plot(yplus,uuplus,'b',label=r'$u^\prime u^\prime$',ls="-")
        ax.plot(yplus,uvplus,'y',label=r'$u^\prime v^\prime$',ls="-")
        ax.plot(yplus,vvplus,'g',label=r'$v^\prime v^\prime$',ls="-")
        ax.plot(yplus,wwplus,'r',label=r'$w^\prime w^\prime$',ls="-")
    ax.minorticks_on()
    ax.set_xscale("symlog",linthresh = 1)
    ax.set_xlabel("$y^+$",fontdict={'size':24})  
    ax.tick_params(axis='x',labelsize=15)
    ax.set_xlim([1,1000])
    x_minor = matplotlib.ticker.LogLocator(base=10.0, subs = np.arange(1.0,10.0))
    ax.xaxis.set_minor_locator(x_minor)
    
    ax.set_ylabel(r'$\xi u^\prime_i u^\prime_j$',\
                  fontdict={'size':24})
    ax.tick_params(axis='y',labelsize=15)

    ax.legend(prop={'size':12},loc='upper right') 
    ax.set_title("Reynolds Stress profile",size=20)

    ax.grid()

    plt.savefig("Reynolds_Stress_profile_wavywall_plusspanwise")
    plt.show()
#%% plot temperature profile
if plt_T : 
    fig, ax = plt.subplots(figsize=[10,8])
    if Compare:
        ax.plot(Tnorm  ,ydelta,  label=r'$T/T_{inf} wavy wall$' ,ls='-')
        ax.plot(Tnorm_f,ydelta_f,label=r'$T/T_{inf} smooth$',ls='--')
    else:
        ax.plot(Tnorm,ydelta,label=r'$T/T_{inf}$')
    ax.minorticks_on()
#    ax.set_xscale("symlog",linthresh = 1)
    ax.set_xlabel(r'$T/T_{inf}$',fontdict={'size':24})  
    ax.tick_params(axis='x',labelsize=15)
    
    ax.set_ylabel(r"$y/\delta$",fontdict={'size':24})
    ax.tick_params(axis='y',labelsize=15)
    ax.set_xlim([1,2])
    ax.set_ylim([0,2])
#    x_minor = matplotlib.ticker.LogLocator(base=10.0, subs = np.arange(1.0,10.0))
#    ax.xaxis.set_minor_locator(x_minor)

    ax.legend(prop={'size':20}) 
    ax.set_title(r"$T/T_{inf}$ profile",size=20)

    ax.grid()

    plt.savefig("temperature_profile_wavywall")
    plt.show()