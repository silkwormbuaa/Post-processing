#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run4.py
@Time    :   2022/11/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Scripts for comparing psd line plot
'''

import os

import sys

source_dir = os.path.realpath(__file__).split('plotlines')[0] 
sys.path.append( source_dir )

import matplotlib.pyplot as     plt

import matplotlib.ticker as     ticker

from   vista.psd         import ProbeData

from   vista.timer       import timer

from   vista.tools       import get_filelist

from   vista.tools       import read_case_parameter


# Option zone
# =============================================================================

output_nr = 5
pure      = False

# =============================================================================


datapath0 = '/home/wencanwu/my_simulation/temp/smooth_wall/probes/probe_x/'
datapath1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/probe_x/'
datapath2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/probes/probe_x/'
datapath3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/probes/probe_x/'
datapath4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/probes/probe_x/'
datapath5 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/probes/probe_x/'

outpath  = '/home/wencanwu/my_simulation/temp/DataPost/'


# smooth wall
with timer("reading smooth wall probe"):
    os.chdir( datapath0 )
    probefile_s  = 'probe_00142.dat'
    Lsep_s       = 13.12627403
    probe_s = ProbeData( probefile_s, withT=False )

if output_nr == 1:
    #1014
    with timer("reading one probe"):
        os.chdir( datapath1 )
        probefile_r = 'probe_00144.dat'
        Lsep_r      = 9.805522
        probe_r = ProbeData( probefile_r, withT=False )
        fig_name = 'psd_line_1014'
        label_r  = r"$D/{\delta}_0=2.0$"

elif output_nr == 2:
    #0926
    with timer("reading one probe"):
        os.chdir( datapath2 )
        probefile_r = 'probe_00145.dat'
        Lsep_r      = 9.676337
        probe_r = ProbeData( probefile_r, withT=False )
        fig_name = 'psd_line_0926'
        label_r  = r"$D/{\delta}_0=1.0$"

elif output_nr == 3:
    #0825
    with timer("reading one probe"):
        os.chdir( datapath3 )
        probefile_r = 'probe_00135.dat'
        Lsep_r      = 11.340638
        probe_r = ProbeData( probefile_r, withT=False )
        fig_name = 'psd_line_0825'
        label_r  = r"$D/{\delta}_0=0.5$"

elif output_nr == 4:
    #0927
    with timer("reading one probe"):
        os.chdir( datapath4 )
        probefile_r = 'probe_00118.dat'
        Lsep_r      = 13.12627
        probe_r = ProbeData( probefile_r, withT=False )
        fig_name = 'psd_line_0927'
        label_r  = r"$D/{\delta}_0=0.25$"

elif output_nr == 5:
    #1221
    with timer("reading one probe"):
        os.chdir( datapath5 )
    #    probefile_r = 'probe_00175.dat'
        probefile_r = 'probe_00192.dat'
        Lsep_r      = 13.266
        probe_r = ProbeData( probefile_r ) 
        fig_name = 'psd_line_1221'   
        label_r  = r"$D/{\delta}_0=0.125$"


""" #1125
with timer("reading one probe"):
    os.chdir( folderpath1 )
#    probefile_r = 'probe_00176.dat'
    probefile_r = 'probe_00189.dat'
    Lsep_r      = 13.12627
    probe_r = ProbeData( probefile_r ) """


with timer("clean data"):
    probe_s.cleandata( starttime=20 )
    probe_r.cleandata( starttime=20 )
    
with timer('psd'):
    
    os.chdir(outpath)
    
    probe_s.psd( 8, 0.3 )
    probe_r.psd( 8, 0.3 )
    
    fig, ax = plt.subplots( figsize=[8,8], constrained_layout=True )
    
    St_Lsep_s = probe_s.St * Lsep_s
    St_Lsep_r = probe_r.St * Lsep_r
    
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

#    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))

    ax.semilogx(St_Lsep_s, 
                probe_s.nd_pprime_fwpsd,
                'red', 
                linewidth=1,
                label='smooth')
    
    ax.semilogx(St_Lsep_r, 
                probe_r.nd_pprime_fwpsd,
                'blue', 
                linewidth=1,
                label=label_r)
    
    ax.set_xlim( [0.01,100] )
    ax.set_ylim( [0.0,1.0] )
    
    ax.set_xlabel( r'$St=f\cdot L_{sep}/U_{\infty}$', 
                  fontdict={'size':24} )
    ax.set_ylabel( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$', 
                  fontdict={'size':24} )

    if not pure:
        ax.legend( prop={'size':20} ) 
    
    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)
        
    if pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
    if pure: fig_name = fig_name + '_pure'
    
    plt.savefig(fig_name)
    
    plt.show()
  
#    print(probe_r.nperseg)
    
#    print(probe_r.df['pprime'])

#    probe_r.write_psd()