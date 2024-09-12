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

output_nr = 6
pure      = True

# =============================================================================


datapath0 = '/media/wencan/Expansion/temp/smooth_isothermal/postprocess/probes/psd_ridge'
datapath1 = '/media/wencan/Expansion/temp/221014/postprocess/probes/psd_ridge'
datapath2 = '/media/wencan/Expansion/temp/220926/postprocess/probes/psd_ridge'
datapath3 = '/media/wencan/Expansion/temp/220825/postprocess/probes/psd_ridge'
datapath4 = '/media/wencan/Expansion/temp/220927/postprocess/probes/psd_ridge'
datapath5 = '/media/wencan/Expansion/temp/221221/postprocess/probes/psd_ridge'
datapath6 = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/probes/psd_ridge'

outpath    = '/media/wencan/Expansion/temp/DataPost/lowRe/psd'


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
    
elif output_nr == 6:
    filelist   = get_filelist( datapath6 )
    
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
            
if output_nr == 0:
    #smooth wall isothermal
    St_Lsep = (np.array(St) * 9.628729).tolist()
    x_sep = [-8.560779, -8.560779]
    y_sep = [0.001, 100]
    x_reatt = [1.06795, 1.06795]
    y_reatt = [0.001, 100]
    x_ppmax = [-7.3333, -7.3333]
    y_ppmax = [0.001, 100]
    fig_name = "psd_smooth_iwall"

elif output_nr == 1:
    #1014
    St_Lsep = (np.array(St) * 9.805522).tolist()
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
    x_sep = [-10.6786859, -10.6786859]
    y_sep = [0.001,100]
    x_reatt = [2.58732, 2.58732]
    y_reatt = [0.001,100]
    x_ppmax = [-8.93096, -8.93096]
    y_ppmax = [0.001, 100]
    fig_name = 'psd_1221'

elif output_nr == 6:
    # smooth adiabatic
    St_Lsep = (np.array(St) * 9.52).tolist()
    x_sep = [-8.42, -8.42]
    y_sep = [0.001,100]
    x_reatt = [1.10, 1.10]
    y_reatt = [0.001,100]
    x_ppmax = [-7.1359675, -7.1359675]
    y_ppmax = [0.001, 100]
    fig_name = "psd_smooth_awall"


# =============================================================================
# plot

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
#plt.show()

