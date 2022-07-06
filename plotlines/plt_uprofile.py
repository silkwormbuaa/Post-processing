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
    y = []
    u = []
    if filename.endswith(".dat"):
        with open(filename) as f:
            for line in f.readlines():
                cleanline = line.strip().split(" ")
                y.append(float(cleanline[0]))
                u.append(float(cleanline[1]))
        fig_label = filename.strip(".dat")    
        print(fig_label)    
        plt.plot(u,y,label=fig_label+"wavy")
lst = os.listdir(".")
lst.sort()
for filename in lst:
    y = []
    u = []
    if filename.endswith("flat.dat"):
        with open(filename) as f:
            for line in f.readlines():
                cleanline = line.strip().split(" ")
                y.append(float(cleanline[0]))
                u.append(float(cleanline[1]))
        fig_label, ext = os.path.splitext(filename)   
        print(fig_label)    
        plt.plot(u,y,label=fig_label,ls="--")
plt.xlabel("<u>",fontdict={'size':24})  
plt.xticks(size = 20)
plt.ylabel("y",fontdict={'size':24})  
plt.yticks(size = 20)  
plt.legend(prop={'size':20})        
plt.grid()
plt.show()
