# -*- coding: utf-8 -*-
'''
@File    :   plt_u.py
@Time    :   2022/08/15 22:59:19
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   plot the velocity profile 
'''
from tkinter import Y
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

#%% Read in data first
FoldPath = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/"
OutPath  = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/DataPost/"
ForceFoldPath = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/forces/forces_1408"

plt_u  = False
plt_RS = False
plt_T  = True

#---Read in averaged flow data
os.chdir(FoldPath)
y_ls    = []
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
#%% plot u profile
os.chdir(OutPath)
if plt_u :
    fig, ax = plt.subplots(figsize=[10,8])
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
    ax.set_title(r"$u^+$ profile wavy wall",size=20)

    ax.grid()

    plt.savefig("velocity_profile_wavywall")
    plt.show()

#%% plot Reynolds Stress profile 
if plt_RS :
    fig, ax = plt.subplots(figsize=[10,8])
    ax.plot(yplus,uuplus,label=r'$u^\prime u^\prime$',ls="--")
    ax.plot(yplus,uvplus,label=r'$u^\prime v^\prime$',ls="--")
    ax.plot(yplus,vvplus,label=r'$v^\prime v^\prime$',ls="--")
    ax.plot(yplus,wwplus,label=r'$w^\prime w^\prime$',ls="--")
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

    ax.legend(prop={'size':20}) 
    ax.set_title("Reynolds Stress profile wavy wall",size=20)

    ax.grid()

    plt.savefig("Reynolds_Stress_profile_wavywall")
    plt.show()
#%% plot temperature profile
if plt_T : 
    fig, ax = plt.subplots(figsize=[10,8])
    ax.plot(ydelta,Tnorm,label=r'$T/T_{inf}$')
    ax.minorticks_on()
#    ax.set_xscale("symlog",linthresh = 1)
    ax.set_xlabel(r"$y/\delta$",fontdict={'size':24})  
    ax.tick_params(axis='x',labelsize=15)
    
    ax.set_ylabel(r'$T/T_{inf}$',fontdict={'size':24})
    ax.tick_params(axis='y',labelsize=15)
    ax.set_xlim([0,2])
#    x_minor = matplotlib.ticker.LogLocator(base=10.0, subs = np.arange(1.0,10.0))
#    ax.xaxis.set_minor_locator(x_minor)

    ax.legend(prop={'size':20}) 
    ax.set_title(r"$T/T_{inf}$ profile wavy wall",size=20)

    ax.grid()

    plt.savefig("temperature_profile_wavywall")
    plt.show()