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

from   plt_tools         import PlotDataframe

import matplotlib.pyplot as     plt

import matplotlib

plt.rcParams.update({'font.size': 18})

OutPath  = "/media/wencanwu/Seagate Expansion Drive/temp/221221/DataPost/"

data0 = '/home/wencanwu/my_simulation/temp/smooth_wall/line_p_prime.dat'

data1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/streamwise_line_1014.dat'

data2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/streamwise_line_0926.dat'

data3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/streamwise_line_0825.dat'

data4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/streamwise_line_0927.dat'

data5 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/streamwise_line_1221.dat'



d1 = PlotDataframe(data1)
d2 = PlotDataframe(data2)
d3 = PlotDataframe(data3)
d4 = PlotDataframe(data4)
d5 = PlotDataframe(data5)
d0 = PlotDataframe(data0)

fig, ax = plt.subplots(figsize=[15,10], constrained_layout=True)

ax.plot( d5.df['(x-x_imp)/δ'], 
         d5.df['<p`>_'],
         'purple',
         label=r'$D/\delta_0=0.125$', 
         ls='-',
#         marker = '^',
#         markevery = 30,
         markersize = 10,
         linewidth=4)
    
ax.plot( d4.df['(x-x_imp)/δ'], 
         d4.df['<p`>_'],
         'red',
         label=r'$D/\delta_0=0.25$', 
         ls     = (0, (10, 3)),
         linewidth=4)

ax.plot( d3.df['(x-x_imp)/δ'], 
         d3.df['<p`>_'],
         'black', 
         label=r'$D/\delta_0=0.5$', 
         ls   = (0, (3, 1, 1, 1, 1, 1)),
         linewidth=4)

ax.plot( d2.df['(x-x_imp)/δ'], 
         d2.df['<p`>_'],
         'blue',  
         label=r'$D/\delta_0=1.0$', 
         ls='-.',
         linewidth=4)

ax.plot( d1.df['(x-x_imp)/δ'], 
         d1.df['<p`>_'],
         'green',
         label=r'$D/\delta_0=2.0$', 
         ls=':',
         linewidth=4)

ax.plot( d0.df['(x-x_imp)/δ'], 
         d0.df['<p`>_'],
         'gray', 
         label=r'$smooth$', 
         ls   ='--',
         linewidth=4)

ax.plot( [-10.6786859,-10.6786859],
         [0.0,0.1],
         'purple',
         linewidth=1)

ax.plot( [-10.6574429,-10.6574429],
         [0.0,0.1],
            'red',
            linewidth=1)

ax.plot( [-9.179693795,-9.179693795],
         [0.0,0.1],
            'black',
            linewidth=1)

ax.plot( [-8.316817364,-8.316817364],
         [0.0,0.1],
            'blue',
            linewidth=1)

ax.plot( [-8.405281852,-8.405281852],
         [0.0,0.1],
            'green',
            linewidth=1)

ax.plot( [-8.56077913,-8.56077913],
         [0.0,0.1],
            'gray',
            linewidth=1)

ax.minorticks_on()

ax.set_xlabel("$(x-x_{imp})/δ_0$",fontdict={'size':40})  
ax.tick_params(axis='x',labelsize=32)

ax.set_ylabel(r"$\sqrt{<p'p'>}/p_{\infty}$",
              fontdict={'size':40})
ax.tick_params(axis='y',labelsize=32)
ax.set_xlim([-20,10])
ax.set_ylim([0.02,0.10])


ax.legend(prop={'size':20}) 


ax.grid()

os.chdir(OutPath)
plt.savefig("pressure_fluctuation_mean")
plt.show()