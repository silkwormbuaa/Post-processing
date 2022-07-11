# -*- coding: utf-8 -*-
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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
fig, ax = plt.subplots(figsize=[10,8])
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
        ax.plot(Tplus,ydelta,label=r'$T^+$',ls="-")


ax.minorticks_on()
#ax.set_xscale("symlog",linthresh = 1)        
ax.set_xlabel(r"$T^+$",fontdict={'size':24})  
ax.tick_params(axis='x',labelsize=15)
#ax.set_xlim([1,1000])
#x_minor = matplotlib.ticker.LogLocator(base=10.0, subs = np.arange(1.0,10.0))
#ax.xaxis.set_minor_locator(x_minor)

ax.set_ylabel(r'$y/\delta$',\
              fontdict={'size':24})
ax.tick_params(axis='y',labelsize=15)
ax.set_xlim([0.5,2.0])
ax.set_ylim([0.0,2.0])
ax.legend(prop={'size':20}) 
ax.set_title(r"$T^+$ profile wavy wall",size=20)   

ax.grid()

plt.savefig("temperature_wavy_wall")
plt.show()