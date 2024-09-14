#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_pressure_fluc.py
@Time    :   2023/09/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   plot pressure(t)
'''


import os
import sys

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

import numpy             as     np
import matplotlib.pyplot as     plt
from   vista.timer       import timer
from   vista.probe       import ProbeData


# Option zone
# =============================================================================

output_nr = 5

# =============================================================================


datapath0 = '/home/wencanwu/my_simulation/temp/smooth_wall/probes/probe_x/'
datapath1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/probe_x/'
datapath2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/probes/probe_x/'
datapath3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/probes/probe_x/'
datapath4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/probes/probe_x/'
datapath5 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/probes/probe_x/'

outpath  = '/home/wencanwu/my_simulation/temp/DataPost/'


with timer("read in probe data "):
    
    if output_nr == 0:
        os.chdir( datapath0 )
        probefile  = 'probe_00142.dat'
    #    probefile = 'probe_00156.dat'
        probe = ProbeData( probefile, withT=False )
        fig_name = 'pressure_0_smooth'
        label  = 'smooth wall'
        
    elif output_nr == 1:
        #1014
        os.chdir( datapath1 )
        probefile = 'probe_00144.dat'
    #    probefile = 'probe_00152.dat'
        probe = ProbeData( probefile, withT=False )
        fig_name = 'pressure_1_1014'
        label  = r"$D/{\delta}_0=2.0$"

    elif output_nr == 2:
        #0926
        os.chdir( datapath2 )
        probefile = 'probe_00145.dat'
    #    probefile = 'probe_00154.dat'
        probe = ProbeData( probefile, withT=False )
        fig_name = 'pressure_2_0926'
        label  = r"$D/{\delta}_0=1.0$" 

    elif output_nr == 3:
        #0825
        os.chdir( datapath3 )
        probefile = 'probe_00135.dat'
    #    probefile = 'probe_00140.dat'
        probe = ProbeData( probefile, withT=False )
        fig_name = 'pressure_3_0825'
        label  = r"$D/{\delta}_0=0.5$"

    elif output_nr == 4:
        #0927
        os.chdir( datapath4 )
        probefile = 'probe_00118.dat'
    #    probefile = 'probe_00125.dat'
        probe = ProbeData( probefile, withT=False )
        fig_name = 'pressure_4_0927'
        label  = r"$D/{\delta}_0=0.25$"
    
    elif output_nr == 5:
        #1221
        os.chdir( datapath5 )
        probefile = 'probe_00175.dat'
    #    probefile = 'probe_00195.dat'
        probe = ProbeData( probefile ) 
        fig_name = 'pressure_5_1221'   
        label  = r"$D/{\delta}_0=0.125$"
    


with timer("clean data"):
    
    probe.cleandata( starttime=20.0 )
    
print( probe.df )

fig, ax = plt.subplots( figsize=[20,5], constrained_layout=True)

time = np.array( probe.df['time'] )
pressure = np.array( probe.df['p'] )/45447.289

ax.plot(time[::20],pressure[::20])

ax.set_ylim( [1.0,2.0])

ax.grid(visible=True)

os.chdir(outpath)
plt.savefig(fig_name+'_separation')
plt.show()