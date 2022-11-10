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

plt.rcParams.update({'font.size': 18})

OutPath  = "/home/wencanwu/my_simulation/temp/221014_lowRe/DataPost/"

data0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/line_p_prime.dat'

data1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/streamwise_line_1014.dat'

data2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/streamwise_line_0926.dat'

data3 = '/home/wencanwu/my_simulation/temp/220825_lowRe/streamwise_line_0825.dat'

data4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/streamwise_line_0927.dat'



d1 = PlotDataframe(data1)
d2 = PlotDataframe(data2)
d3 = PlotDataframe(data3)
d4 = PlotDataframe(data4)
d0 = PlotDataframe(data0)

fig, ax = plt.subplots(figsize=[10,8.5], constrained_layout=True)
    
ax.plot( d0.df['(x-x_imp)/δ'], 
         d0.df['<p`>_'],
         'gray', 
         label=r'$smooth$', 
         ls   ='--')

ax.plot( d1.df['(x-x_imp)/δ'], 
         d1.df['<p`>_'],
         'green',
         label=r'$D=2.0\delta_0$', 
         ls='-',
         marker = 's',
         markevery = 28,
         markersize = 10)

ax.plot( d2.df['(x-x_imp)/δ'], 
         d2.df['<p`>_'],
         'blue',  
         label=r'$D=1.0\delta_0$', 
         ls='-',
         marker = 'o',
         markevery = 29,
         markersize = 10)

ax.plot( d3.df['(x-x_imp)/δ'], 
         d3.df['<p`>_'],
         'black', 
         label=r'$D=0.5\delta_0$', 
         ls='-',
         marker = 'v',
         markevery = 30,
         markersize = 10)

ax.plot( d4.df['(x-x_imp)/δ'], 
         d4.df['<p`>_'],
         'red',
         label=r'$D=0.25\delta_0$', 
         ls='-',
         marker = 'D',
         markevery = 30,
         markersize = 10)

ax.minorticks_on()

ax.set_xlabel("$(x-x_{imp})/δ_0$",fontdict={'size':40})  
ax.tick_params(axis='x',labelsize=32)

ax.set_ylabel(r"$\sqrt{<p'p'>}/p_{\infty}$",
              fontdict={'size':40})
ax.tick_params(axis='y',labelsize=32)
ax.set_xlim([-20,10])


#ax.legend(prop={'size':28}) 


ax.grid()

os.chdir(OutPath)
plt.savefig("pressure_fluctuation_mean_marker")
plt.show()