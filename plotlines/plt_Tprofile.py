# -*- coding: utf-8 -*-
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
import os
"""
Before run, modify the variable for x axis 
"""

#%% plot streamwise velocity profile
OutFile  = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/DataPost/"
#OutFile  = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/DataPost/"
os.chdir(OutFile)
lst = os.listdir(OutFile)
lst.sort()
plt.figure(figsize=[12,10])
for filename in lst:
    yplus  = []
    yplus2 = []
    y      = []
    T  = []
    if filename.endswith("wavy.dat"):
        with open(filename) as f:
            for line in f.readlines():
                cleanline = line.strip().split(" ")
                yplus.append(float(cleanline[0]))
                T.append(float(cleanline[1]))
                yplus2.append(float(cleanline[2]))
                y.append(float(cleanline[3]))                
        fig_label = filename.strip(".dat")    
        print(fig_label)    
        plt.plot(yplus2,T,label=fig_label)
lst = os.listdir(OutFile)
lst.sort()
for filename in lst:
    y = []
    T = []
    yplus = []
    yplus2= []    
    if filename.endswith("flat.dat"):
        with open(filename) as f:
            for line in f.readlines():
                cleanline = line.strip().split(" ")
                yplus.append(float(cleanline[0]))
                T.append(float(cleanline[1]))
                yplus2.append(float(cleanline[2]))
                y.append(float(cleanline[3]))
        fig_label, ext = os.path.splitext(filename)   
        print(fig_label)    
        plt.plot(yplus2,T,label=fig_label,ls="--")       
plt.xscale("symlog")        
plt.ylabel("T",fontdict={'size':24})  
plt.xticks(size = 20)
plt.xlabel("$y^+$",fontdict={'size':24})  
plt.yticks(size = 20)  
plt.legend(prop={'size':20}) 
plt.title("Temperature profile")      
#plt.xlim(0,1000)
plt.grid()
plt.show()
