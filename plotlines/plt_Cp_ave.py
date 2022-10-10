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

FoldPath = '/home/wencanwu/my_simulation/temp/220927_lowRe/linedata'

os.chdir(FoldPath)

delta = 5.2
x_imp = 50.4
x_hd  = []
Cf_hd = []
Pf_hd = []
with open('Cf_points_0825.dat') as f:
    for line in f.readlines():
        cleanl = line.strip().split()
        x_hd .append((float(cleanl[0])-x_imp)/delta)
        Cf_hd.append(float(cleanl[1])*1000)
        Pf_hd.append(float(cleanl[2]))


x_flat  = []
Cf_flat = []
Pf_flat = []
with open('Cf_flat_new.dat') as f:
    for line in f.readlines():
        cleanl = line.strip().split()
        x_flat. append((float(cleanl[0])-x_imp)/delta)
        Cf_flat.append(float(cleanl[1])*1000)
        Pf_flat.append(float(cleanl[2]))

x_1d  = []
Cf_1d = []
Pf_1d = []
with open('Cf_points_0926.dat') as f:
    for line in f.readlines():
        cleanl = line.strip().split()
        x_1d. append((float(cleanl[0])-x_imp)/delta)
        Cf_1d.append(float(cleanl[1])*1000)
        Pf_1d.append(float(cleanl[2]))

x_qd  = []
Cf_qd = []
Pf_qd = []
with open('Cf_points_0927.dat') as f:
    for line in f.readlines():
        cleanl = line.strip().split()
        x_qd. append((float(cleanl[0])-x_imp)/delta)
        Cf_qd.append(float(cleanl[1])*1000)
        Pf_qd.append(float(cleanl[2]))
'''
fig, ax = plt.subplots(figsize=[10,8])
ax.plot(x_flat,Cf_flat,'gray',label=r'$smooth$',ls='--')
ax.plot(x_crest,Cf_crest,'b',label=r'$z=0.0$',ls='-')
ax.plot(x_valley,Cf_valley,'r',label=r'$z=1.3$',ls='-')
ax.plot(x_point,Cf_point,'black',label=r'$averaged$', marker="s")
ax.minorticks_on()

ax.set_xlabel("$(x-x_{imp})/\delta_0$",fontdict={'size':24})
ax.tick_params(axis='x',labelsize=15)
ax.set_ylabel("$C_fx10^3$",fontdict={'size':24})
ax.tick_params(axis='y',labelsize=15)
ax.set_xlim([-20.0,15.0])
ax.set_ylim([-1.0,4.0])

ax.legend(prop={'size':15})
ax.grid()

plt.savefig("Cf.png")
plt.show()

'''
fig, ax = plt.subplots(figsize=[10,8])
ax.plot(x_flat,Pf_flat,'gray',label=r'$smooth$',ls='--')
ax.plot(x_1d,Pf_1d,'b',label=r'$D=1.0\delta_0$',marker = 's')
ax.plot(x_hd,Pf_hd,'black',label=r'$D=0.5\delta_0$',marker = 's')
ax.plot(x_qd,Pf_qd,'r',label=r'$D=0.25\delta_0$',marker = 's')
ax.minorticks_on()

ax.set_xlabel("$(x-x_{imp})/\delta_0$",fontdict={'size':24})
ax.tick_params(axis='x',labelsize=15)
ax.set_ylabel("$p_w/p_{\infty}$",fontdict={'size':24})
ax.tick_params(axis='y',labelsize=15)
ax.set_xlim([-20.0,15.0])
#ax.set_ylim([0.3,0.85])

ax.legend(prop={'size':15})
ax.grid()

plt.savefig("C_p.png")
plt.show()
