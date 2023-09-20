# -*- coding: utf-8 -*-
'''
@File    :   plt_psd_pure.py
@Time    :   2023/02/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plotting psd contours with/without labels
'''

import os

import sys

import numpy             as     np

import matplotlib        as     mpl

import matplotlib.pyplot as     plt

from   matplotlib        import cm

from   matplotlib.colors import ListedColormap,LinearSegmentedColormap

from   plt_tools         import *

source_dir = os.path.realpath(__file__).split('plotlines')[0] 
sys.path.append( source_dir )

from   vista.tools       import get_filelist

from   vista.timer       import timer

plt.rcParams.update({'font.size': 25})

# Option zone
# =============================================================================

output_nr = 5
pure      = True

# =============================================================================


datapath0 = '/home/wencanwu/my_simulation/temp/smooth_wall/probes/psd_x'
datapath1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/psd_x'
datapath2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/probes/psd_x'
datapath3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/probes/psd_x'
datapath4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/probes/psd_x'
datapath5 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/probes/psd_x'

outpath    = '/home/wencanwu/my_simulation/temp/DataPost'


if output_nr == 0:
    filelist   = get_filelist( datapath0 )

elif output_nr == 1:
    filelist   = get_filelist( datapath1 )

elif output_nr == 2:
    filelist   = get_filelist( datapath2 )

elif output_nr == 3:
    filelist   = get_filelist( datapath3 )

elif output_nr == 4:
    filelist   = get_filelist( datapath4 )

elif output_nr == 5:
    filelist   = get_filelist( datapath5 )
    
else: raise ValueError("case nr is wrong!")
    
x = None
St = None
nd_fwpsd = None

with timer("import psd data"):
    for i, psdfile in enumerate( filelist ):
        
        print(i, psdfile[-9:])
        with timer("read one psd file"):
            x_temp, St_temp, nd_fwpsd_temp = read_psd( psdfile )
            print(len(St_temp))
        
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
            
if output_nr == 0:
    #smooth wall
    St_Lsep = (np.array(St) * 9.628729).tolist()
    print(St_Lsep[0][:10])
    print(St_Lsep[0][-10:])
    x_sep = [-8.560779, -8.560779]
    y_sep = [0.001, 100]
    x_reatt = [1.06795, 1.06795]
    y_reatt = [0.001, 100]
    x_ppmax = [-7.3333, -7.3333]
    y_ppmax = [0.001, 100]
    fig_name = "psd_smooth"

elif output_nr == 1:
    #1014
    St_Lsep = (np.array(St) * 9.805522).tolist()
    print(St_Lsep[0][:10])
    x_sep = [-8.4052818, -8.4052818]
    y_sep = [0.001,100]
    x_reatt = [1.4002403, 1.4002403]
    y_reatt = [0.001,100]
    x_ppmax = [-7.712139, -7.712139]
    y_ppmax = [0.001, 100]
    fig_name = "psd_1014"

elif output_nr == 2:
    #0926
    St_Lsep = (np.array(St) * 9.676337).tolist()
    print(St_Lsep[0][:10])
    x_sep = [-8.316817, -8.316817]
    y_sep = [0.001,100]
    x_reatt = [1.3595199, 1.3595199]
    y_reatt = [0.001,100]
    x_ppmax = [-7.49053, -7.49053]
    y_ppmax = [0.001, 100]
    fig_name = "psd_0926"

elif output_nr == 3:
     #0825
    St_Lsep = (np.array(St) * 11.340638).tolist()
    print(St_Lsep[0][:10])
    x_sep = [-9.1796937, -9.1796937]
    y_sep = [0.001,100]
    x_reatt = [2.1609445, 2.1609445]
    y_reatt = [0.001,100]
    x_ppmax = [-8.70936, -8.70936]
    y_ppmax = [0.001, 100]
    fig_name = "psd_0825"

elif output_nr == 4:
    #0927
    St_Lsep = (np.array(St) * 13.12627).tolist()
    print(St_Lsep[0][:10])
    print(St_Lsep[0][-10:])
    x_sep = [-10.65744, -10.65744]
    y_sep = [0.001,100]
    x_reatt = [2.468831, 2.468831]
    y_reatt = [0.001,100]
    x_ppmax = [-10.06114, -10.06114]
    y_ppmax = [0.001, 100]
    fig_name = "psd_0927"

elif output_nr == 5:
    ##1221
    St_Lsep = (np.array(St) * 13.266).tolist()
    print(St_Lsep[0][:10])
    print(St_Lsep[0][-10:])
    x_sep = [-10.6786859, -10.6786859]
    y_sep = [0.001,100]
    x_reatt = [2.58732, 2.58732]
    y_reatt = [0.001,100]
    x_ppmax = [-8.93096, -8.93096]
    y_ppmax = [0.001, 100]
    fig_name = 'psd_1221'
    
    #1125
#    St_Lsep = (np.array(St) * 13.1286891).tolist()
#    print(St_Lsep[0][:10])
#    print(St_Lsep[0][-10:])
#    x_sep = [-10.55859, -10.55859]
#    y_sep = [0.001,100]
#    x_reatt = [2.570097, 2.570097]
#    y_reatt = [0.001,100]
#    fig_name = "psd_1125"



fig, ax = plt.subplots( figsize=[15,6],
                        constrained_layout=True)

fwpsd = ax.pcolormesh(x,St_Lsep,nd_fwpsd,
                      shading='gouraud',
                      cmap='Greys',
                      vmin=0, vmax=0.3)

# -- using contourf instead of pcolormesh

#fwpsd = ax.contourf(x,St_Lsep,nd_fwpsd,np.linspace(0,0.3,num=31),
#                      cmap='Greys',vmin=0, vmax=0.3)

if not pure:
    fig.colorbar( fwpsd ,orientation='horizontal')

ax.plot(x_sep,y_sep,'r',linewidth=1.0)
ax.plot(x_reatt,y_reatt,'r',linewidth=1.0)
ax.plot(x_ppmax,y_ppmax,'blue',linewidth=1.0)

ax.minorticks_on()
ax.tick_params(which='major',
                axis='both',
                direction='out',
                length=20,
                width=2.0)
ax.tick_params(which='minor',
                axis='both', 
                direction='out',
                length=10,
                width=1.5)

ax.set_yscale('log')
ax.set_ylim( [0.01,100] )

ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

if not pure:
    ax.set_xlabel( r'$(x-x_{imp})/\delta_0$', fontdict={'size':16})
    ax.set_ylabel( r'$St=f\cdot L_{sep}/U_{\infty}$', fontdict={'size':16})

    ax.set_title( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$', size=18)

os.chdir( outpath )

if pure: fig_name = fig_name + '_pure_ppmax'
plt.savefig( fig_name )
#plt.show()

