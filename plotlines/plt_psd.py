#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_psd.py
@Time    :   2022/11/14 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plotting PSD contours
'''

import os

import sys

import numpy             as     np

import matplotlib        as     mpl

import matplotlib.pyplot as     plt

from   matplotlib        import cm

from   matplotlib.colors import ListedColormap,LinearSegmentedColormap

from   plt_tools         import *

sys.path.append("..") 

from   vista_tools       import *

from   timer             import timer

plt.rcParams.update({'font.size': 16})

folderpath = '/home/wencanwu/my_simulation/temp/220927_lowRe/probes/psd_x'
outpath    = '/home/wencanwu/my_simulation/temp/220927_lowRe/probes'

filelist   = get_filelist(folderpath)

os.chdir( folderpath )

x = None
St = None
nd_fwpsd = None

with timer("import psd data"):
    for i, psdfile in enumerate( filelist ):
        
        print(i)
        with timer("read one psd file"):
            x_temp, St_temp, nd_fwpsd_temp = read_psd( psdfile )
        
        if i == 0:
            
            x = list()
            St = list()
            nd_fwpsd = list()
            x.append(x_temp)
            St.append(St_temp)
            nd_fwpsd.append(nd_fwpsd_temp)
        
        else:
            
            x.append(x_temp)
            St.append(St_temp)
            nd_fwpsd.append(nd_fwpsd_temp)

#1014
#   St_Lsep = (np.array(St) * 9.805522).tolist()
#   print(St_Lsep[0][:10])
#   x_sep = [-8.4052818, -8.4052818]
#   y_sep = [0.001,100]
#   x_reatt = [1.4002403, 1.4002403]
#   y_reatt = [0.001,100]

#    #0926
#    St_Lsep = (np.array(St) * 9.676337).tolist()
#    print(St_Lsep[0][:10])
#    x_sep = [-8.316817, -8.316817]
#    y_sep = [0.001,100]
#    x_reatt = [1.3595199, 1.3595199]
#    y_reatt = [0.001,100]
#    #0825
#   St_Lsep = (np.array(St) * 11.340638).tolist()
#   print(St_Lsep[0][:10])
#   x_sep = [-9.1796937, -9.1796937]
#   y_sep = [0.001,100]
#   x_reatt = [2.1609445, 2.1609445]
#   y_reatt = [0.001,100]
##0927
St_Lsep = (np.array(St) * 13.12627).tolist()
print(St_Lsep[0][:10])
print(St_Lsep[0][-10:])
x_sep = [-10.65744, -10.65744]
y_sep = [0.001,100]
x_reatt = [2.468831, 2.468831]
y_reatt = [0.001,100]


top = cm.get_cmap('Blues_r', 128) # r means reversed version
bottom = cm.get_cmap('Oranges', 128)  # combine it all
newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
orange_blue = ListedColormap(newcolors, name='OrangeBlue')

fig, ax = plt.subplots( figsize=[9,4],
                       constrained_layout=True)

fwpsd = ax.pcolormesh(x,St_Lsep,nd_fwpsd,
                      shading='gouraud',
                      cmap='Greys',
                      vmin=0, vmax=0.3)
fig.colorbar( fwpsd )

ax.plot(x_sep,y_sep,'r',linewidth=0.3)
ax.plot(x_reatt,y_reatt,'r',linewidth=0.3)

ax.minorticks_on()
ax.set_xlabel( r'$(x-x_{imp})/\delta_0$', fontdict={'size':16})
ax.set_ylabel( r'$St=f\cdot L_{sep}/U_{\infty}$', fontdict={'size':16})
ax.set_yscale('log')
ax.set_ylim( [0.001,100] )

ax.set_title( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$',
             size=18)

os.chdir(os.pardir)
plt.savefig('psd0927_Lsep_gray.png')
#plt.show()

