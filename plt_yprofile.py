# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os

#%% plot streamwise velocity profile
lst = os.listdir(".")
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
        plt.plot(u,y,label=filename)
      
plt.legend()        
plt.grid()
plt.show()
