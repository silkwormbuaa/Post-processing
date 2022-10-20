#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   pli_p_prime.py
@Time    :   2022/10/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os

from   plt_tools         import *

import matplotlib.pyplot as     plt

import matplotlib

OutPath  = "/home/wencanwu/my_simulation/temp/220927_lowRe/DataPost/"

data0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/line_p_prime.dat'

data1 = '/home/wencanwu/my_simulation/temp/220926_lowRe/line_crest.dat'

data2 = '/home/wencanwu/my_simulation/temp/220825_lowRe/line_crest.dat'

data3 = '/home/wencanwu/my_simulation/temp/220927_lowRe/line_crest.dat'

d1 = PlotDataframe(data1)
d2 = PlotDataframe(data2)
d3 = PlotDataframe(data3)
d0 = PlotDataframe(data0)

fig, ax = plt.subplots(figsize=[16,8])
    
ax.plot(d0.df['(x-x_imp)/δ'], d0.df['<p`>_'],
            'gray', label=r'$smooth$', ls='--')

ax.plot(d1.df['(x-x_imp)/δ'], d1.df['<p`>_'],
        'blue',  label=r'$D=1.0\delta_0$', ls='-')

ax.plot(d2.df['(x-x_imp)/δ'], d2.df['<p`>_'],
        'black', label=r'$D=0.5\delta_0$', ls='-')

ax.plot(d3.df['(x-x_imp)/δ'], d3.df['<p`>_'],
        'red',   label=r'$D=0.25\delta_0$', ls='-')

ax.minorticks_on()

ax.set_xlabel("$(x-x_{imp})/δ$",fontdict={'size':24})  
ax.tick_params(axis='x',labelsize=15)

ax.set_ylabel('$p\'/p_{\infty}$',fontdict={'size':24})
ax.tick_params(axis='y',labelsize=15)
ax.set_xlim([-20,10])


ax.legend(prop={'size':20}) 
ax.set_title(r"$u^+$ profile",size=20)

ax.grid()

os.chdir(OutPath)
plt.savefig("pressure_fluctuation_crest")
plt.show()