# -*- coding: utf-8 -*-
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
import os

#%% plot streamwise velocity profile
OutFile  = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/DataPost/"
#OutFile  = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/DataPost/"
os.chdir(OutFile)
lst = os.listdir(OutFile)
lst.sort()
plt.figure(figsize=[12,10])
for filename in lst:
    y  = []
    uu = []
    uv = []
    vv = []
    ww = []
    if filename.endswith("wavy.dat"):
        with open(filename) as f:
            for line in f.readlines():
                cleanline = line.strip().split(" ")
                y.append(float(cleanline[0]))
                uu.append(float(cleanline[1]))
                uv.append(float(cleanline[2]))
                vv.append(float(cleanline[3]))
                ww.append(float(cleanline[4]))
        fig_label = filename.strip(".dat")    
        print(fig_label)    
        plt.plot(y,uu,label=r'$u^\prime u^\prime$')
        plt.plot(y,uv,label=r'$u^\prime v^\prime$')
        plt.plot(y,vv,label=r'$v^\prime v^\prime$')
        plt.plot(y,ww,label=r'$w^\prime w^\prime$')
lst = os.listdir(OutFile)
lst.sort()
for filename in lst:
    y = []
    uu = []
    uv = []
    vv = []
    ww = []
    if filename.endswith("flat.dat"):
        with open(filename) as f:
            for line in f.readlines():
                cleanline = line.strip().split(" ")
                y.append(float(cleanline[0]))
                uu.append(float(cleanline[1]))
                uv.append(float(cleanline[2]))
                vv.append(float(cleanline[3]))
                ww.append(float(cleanline[4]))
        fig_label, ext = os.path.splitext(filename)   
        print(fig_label)    
        plt.plot(y,uu,label=r'$u^\prime u^\prime$',ls="--")
        plt.plot(y,uv,label=r'$u^\prime v^\prime$',ls="--")
        plt.plot(y,vv,label=r'$v^\prime v^\prime$',ls="--")
        plt.plot(y,ww,label=r'$w^\prime w^\prime$',ls="--")        
plt.xscale("symlog")        
plt.ylabel(r"$\sqrt{\langle u^{\prime}_i u^{\prime}_j\rangle^{+}}$",fontdict={'size':24})  
plt.xticks(size = 20)
plt.xlabel("$y^+$",fontdict={'size':24})  
plt.yticks(size = 20)  
plt.legend(prop={'size':20}) 
plt.title(fig_label)      
plt.xlim(0,1000)
plt.grid()
plt.show()
