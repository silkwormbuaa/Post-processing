# -*- coding: utf-8 -*-
'''
@File    :   plt_Cf.py
@Time    :   2022/08/30 18:27:00
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

FoldPath = '/home/wencanwu/my_simulation/temp/220825_lowRe/linedata'

os.chdir(FoldPath)
delta = 5.2
x_flat  = []
Cf_flat = []
Pf_flat = []
with open('x_cf_STBLI_Wencan.dat') as f:
    for line in f.readlines():
        cleanl = line.strip().split()
        x_flat. append(float(cleanl[0])/delta)
        Cf_flat.append(float(cleanl[1])*1000)
#        Pf_flat.append(float(cleanl[2]))

x_crest  = []
Cf_crest = []
Pf_crest = []
with open('Cf_crest.dat') as f:
    for line in f.readlines():
        cleanl = line.strip().split()
        x_crest. append(float(cleanl[0])/delta)
        Cf_crest.append(float(cleanl[1])*1000)
        Pf_crest.append(float(cleanl[2]))

x_valley  = []
Cf_valley = []
Pf_valley = []
with open('Cf_valley.dat') as f:
    for line in f.readlines():
        cleanl = line.strip().split()
        x_valley. append(float(cleanl[0])/delta)
        Cf_valley.append(float(cleanl[1])*1000)
        Pf_valley.append(float(cleanl[2]))

fig, ax = plt.subplots(figsize=[10,8])
ax.plot(x_flat,Cf_flat,'gray',label=r'$smooth$',ls='--')
ax.plot(x_crest,Cf_crest,'b',label=r'$z=0.0$',ls='-')
ax.plot(x_valley,Cf_valley,'r',label=r'$z=1.3$',ls='-')
ax.minorticks_on()

ax.set_xlabel("$x/\delta_0$",fontdict={'size':24})
ax.tick_params(axis='x',labelsize=15)
ax.set_ylabel("$C_fx10^3$",fontdict={'size':24})
ax.tick_params(axis='y',labelsize=15)
ax.set_xlim([-5.0,20.0])

ax.legend(prop={'size':15})
ax.grid()

plt.savefig("Cf.png")
plt.show()

'''
fig, ax = plt.subplots(figsize=[10,8])
ax.plot(x_flat,Pf_flat,'gray',label=r'$smooth$',ls='--')
ax.plot(x_crest,Pf_crest,'b',label=r'$z=0.0$',ls='-')
ax.plot(x_valley,Pf_valley,'r',label=r'$z=1.3$',ls='-')
ax.minorticks_on()

ax.set_xlabel("$x/\delta_0$",fontdict={'size':24})
ax.tick_params(axis='x',labelsize=15)
ax.set_ylabel("$p_w/p_{\infty}$",fontdict={'size':24})
ax.tick_params(axis='y',labelsize=15)
ax.set_xlim([-5.0,20.0])

ax.legend(prop={'size':15})
ax.grid()

plt.savefig("Pf.png")
plt.show()
'''