#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_psd_map.py
@Time    :   2024/03/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.timer       import timer
from   vista.tools       import get_filelist

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = "Times New Roman"

# Option zone
# =============================================================================

output_nr = 1              # 1,2,3
pure      = False

# =============================================================================


datapath1 = '/media/wencanwu/Seagate Expansion Drive1/temp/240211/postprocess/probes/psd_ridge'
datapath2 = '/media/wencanwu/Seagate Expansion Drive1/temp/220927/postprocess/probes/psd_ridge'
datapath3 = '/media/wencanwu/Seagate Expansion Drive1/temp/240210/postprocess/probes/psd_ridge'

outpath    = '/media/wencanwu/Seagate Expansion Drive1/temp/DataPost/lowRe_ridge_height/psd'


if output_nr == 1:
    filelist   = get_filelist( datapath1 )

elif output_nr == 2:
    filelist   = get_filelist( datapath2 )

elif output_nr == 3:
    filelist   = get_filelist( datapath3 )
    
else: raise ValueError("case nr is wrong!")
    
x = None
St = None
pmpsd = None

with timer("import psd data"):
    for i, psdfile in enumerate( filelist ):
        
        print(f"Imported {i:05d} probes: {psdfile[-9:]}")
        with timer("read one psd file"):
            
            probe = ProbeData()
            probe.read_psd( psdfile )
            st_temp    = np.array( probe.psd_df['freq'] )*5.2/507.0
            x_temp     = [(probe.xyz[0]-50.4)/5.2]*len(st_temp)
            pmpsd_temp = np.array( probe.psd_df['pmpsd_p_fluc'] )

        if i == 0:
            
            x = list()
            St = list()
            pmpsd = list()
            x.append(x_temp)
            St.append(st_temp)
            pmpsd.append(pmpsd_temp)
        
        else:
            
            x.append(x_temp)
            St.append(st_temp)
            pmpsd.append(pmpsd_temp)

if output_nr == 1:
    #240211
    St_Lsep = (np.array(St) * 11.44).tolist()
    x_sep = [-9.510, -9.510]
    y_sep = [0.001,100]
    x_reatt = [1.93, 1.93]
    y_reatt = [0.001,100]
    x_ppmax = [-8.7426007, -8.7426007]
    y_ppmax = [0.001, 100]
    fig_name = "psd_240211"

elif output_nr == 2:
    #0927
    St_Lsep = (np.array(St) * 13.12627).tolist()
    x_sep = [-10.65744, -10.65744]
    y_sep = [0.001,100]
    x_reatt = [2.468831, 2.468831]
    y_reatt = [0.001,100]
    x_ppmax = [-10.06114, -10.06114]
    y_ppmax = [0.001, 100]
    fig_name = "psd_220927"

elif output_nr == 3:
    ##240210
    St_Lsep = (np.array(St) * 17.1).tolist()
    x_sep = [-13.85, -13.85]
    y_sep = [0.001,100]
    x_reatt = [3.25, 3.25]
    y_reatt = [0.001,100]
    x_ppmax = [-11.5126578, -11.5126578]
    y_ppmax = [0.001, 100]
    fig_name = 'psd_240210'
    
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

fwpsd = ax.pcolormesh(x,St_Lsep,pmpsd,
                      shading='gouraud',
                      cmap='Greys',
                      vmin=0, vmax=0.3)

# -- using contourf instead of pcolormesh

#fwpsd = ax.contourf(x,St_Lsep,pmpsd,np.linspace(0,0.3,num=31),
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