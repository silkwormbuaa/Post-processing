#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_p_prime.py
@Time    :   2022/10/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os

from   plt_tools         import PlotDataframe

import matplotlib.pyplot as     plt

import matplotlib.ticker as     ticker

plt.rcParams.update({'font.size': 18})

OutPath  = "/home/wencanwu/my_simulation/temp/DataPost/"

data0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/line_p_prime.dat'

data1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/streamwise_line_1014.dat'

data2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/streamwise_line_0926.dat'

data3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/streamwise_line_0825.dat'

data4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/streamwise_line_0927.dat'

data5 = '/home/wencanwu/my_simulation/temp/221125_lowRe/streamwise_line_1125.dat'



d1 = PlotDataframe(data1)
d2 = PlotDataframe(data2)
d3 = PlotDataframe(data3)
d4 = PlotDataframe(data4)
d5 = PlotDataframe(data5)
d0 = PlotDataframe(data0)

fig, ax = plt.subplots(figsize=[8,8], constrained_layout=True)
    

ax.plot( d2.df['(x-x_imp)/δ'], 
         d2.df['<p`>_'],
         'blue',  
         label=r'$D/\delta_0=1.0$', 
         ls='-.',
         linewidth=4)

ax.plot( d3.df['(x-x_imp)/δ'], 
         d3.df['<p`>_'],
         'black', 
         label=r'$D/\delta_0=0.5$', 
         ls   = (0, (3, 1, 1, 1, 1, 1)),
         linewidth=4)

ax.plot( d4.df['(x-x_imp)/δ'], 
         d4.df['<p`>_'],
         'red',
         label=r'$D/\delta_0=0.25$', 
         ls     = (0, (10, 3)),
         linewidth=4)

#ax.plot( d5.df['(x-x_imp)/δ'], 
#         d5.df['<p`>_'],
#         'purple',
#         label=r'$D/\delta_0=0.125$', 
#         ls='-',
#         marker = '^',
#         markevery = 30,
#         markersize = 10)

ax.plot( d0.df['(x-x_imp)/δ'], 
         d0.df['<p`>_'],
         'gray', 
         label=r'$smooth$', 
         ls   ='--',
         linewidth=4)

ax.plot( d1.df['(x-x_imp)/δ'], 
         d1.df['<p`>_'],
         'green',
         label=r'$D/\delta_0=2.0$', 
         ls=':',
         linewidth=4)

ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))

ax.minorticks_on()
ax.tick_params(which='major',
                axis='both',
                direction='in',
                length=15,
                width=2)
ax.tick_params(which='minor',
                axis='both', 
                direction='in',
                length=10,
                width=1)

ax.set_xlabel("$(x-x_{imp})/δ_0$",fontdict={'size':40})  
ax.tick_params(axis='x',labelsize=32)

ax.set_ylabel(r"$\sqrt{<p'p'>}/p_{\infty}$",
              fontdict={'size':40})
ax.tick_params(axis='y',labelsize=32)

ax.set_xlim([-20,10])
ax.set_ylim([0.02,0.1])

#ax.legend(prop={'size':20}) 

# Adjust the spacing around the plot to remove the white margin
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)

os.chdir(OutPath)
plt.savefig("pressure_fluctuation_pure")
plt.show()