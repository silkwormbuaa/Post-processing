# -*- coding: utf-8 -*-
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
import os
"""
Before run, modify the variable for x axis 
"""

#%% plot streamwise velocity profile
#OutFile  = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/DataPost/"
OutFile  = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/DataPost/"
os.chdir(OutFile)
lst = os.listdir(OutFile)
lst.sort()
fig, ax = plt.subplots()
plt.figure(figsize=[10,8])
for filename in lst:
    xlst   = []
    ylst   = []
    deltay = []
    yplus  = []
    yplus2 = []
    ydelta = []
    ulst   = []
    uplus  = []
    uunor  = []
    vvnor  = []
    wwnor  = []
    uvnor  = []
    Tlst   = []
    Tplus  = []
    if filename.endswith(".dat"):
        with open(filename) as f:
            for line in f.readlines():
                cleanline = line.strip().split(" ")
                xlst.append(float(  cleanline[0]))
                ylst.append(float(  cleanline[1]))  
                deltay.append(float(cleanline[2]))
                yplus.append(float( cleanline[3]))
                yplus2.append(float(cleanline[4])) 
                ydelta.append(float(cleanline[5]))
                ulst.append(float(  cleanline[6]))  
                uplus.append(float( cleanline[7])) 
                uunor.append(float( cleanline[8])) 
                vvnor.append(float( cleanline[9])) 
                wwnor.append(float( cleanline[10])) 
                uvnor.append(float( cleanline[11])) 
                Tlst.append(float(  cleanline[12]))  
                Tplus.append(float( cleanline[13])) 
        fig_label = filename.strip(".dat")    
        print(fig_label)    
        ax.plot(yplus,uunor,label=r'$u^\prime u^\prime$',ls="--")
        ax.plot(yplus,uvnor,label=r'$u^\prime v^\prime$',ls="--")
        ax.plot(yplus,vvnor,label=r'$v^\prime v^\prime$',ls="--")
        ax.plot(yplus,wwnor,label=r'$w^\prime w^\prime$',ls="--")
        
plt.xscale("symlog")        
plt.ylabel(r"$\sqrt{\langle u^{\prime}_i u^{\prime}_j\rangle^{+}}$",\
          fontdict={'size':24})  
plt.xticks(size = 20)
plt.xlabel("$y^+$",fontdict={'size':24})  
plt.yticks(size = 20)  
#plt.ylim([-2,8])
plt.minorticks_on()
#plt.tick_params(direction="out")
plt.legend(prop={'size':20}) 
plt.title("Reynolds Stress profile_Flat Plate")      
plt.xlim(1,1000)
plt.grid()
plt.show()
