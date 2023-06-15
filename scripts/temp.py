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

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

import matplotlib.pyplot as     plt

import matplotlib.ticker as     ticker

from   vista.psd         import ProbeData

from   vista.timer       import timer

from   vista.tools       import get_filelist

from   vista.tools       import read_case_parameter


folderpath0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/probes'
folderpath1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/probe_x'
#folderpath1 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/probes/probe_x'

#outpath0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/probes/psd_x_new'
#outpath1 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/probes/psd_x'
#outpath  = '/home/wencanwu/my_simulation/temp/DataPost/'
outpath1 = '/home/wencanwu/my_simulation/temp/221014_lowRe/probes/psd_x'


Pure = False # if output pure figure for latex
#%%
"""
os.chdir(folderpath1)

filelist = get_filelist( folderpath1 )

n_probe = len(filelist)

print( "we have got %5d probes data"%( n_probe ) )

with timer('PSD for one case'):
    
    os.chdir( outpath1 )
    
    for i, probefile in enumerate(filelist):
        
        with timer('%5.3f%% Get psd for %s probe'%(float(i)/n_probe*100,probefile[-9:-4])):
            
            probe = ProbeData( probefile )
            
            probe.cleandata( starttime=20 )
            
            probe.psd( 8, 0.5 )

            probe.write_psd()

"""
#%% print single line

# smooth wall
with timer("reading smooth wall probe"):
    os.chdir( folderpath0 )
    probefile_s  = 'probe_sep.dat'
    Lsep_s       = 13.12627403
    probe_s = ProbeData( probefile_s )
"""
#1014
with timer("reading one probe"):
    os.chdir( folderpath1 )
    probefile1 = 'probe_00144.dat'
    Lsep1      = 9.805522
    probe1 = ProbeData( probefile1 )
"""
"""    

#0926
with timer("reading one probe"):
    os.chdir( folderpath1 )
    probefile1 = 'probe_00145.dat'
    Lsep1      = 9.676337
    probe1 = ProbeData( probefile1 )

#0825
with timer("reading one probe"):
    os.chdir( folderpath1 )
    probefile1 = 'probe_00135.dat'
    Lsep1      = 11.340638
    probe1 = ProbeData( probefile1 )
#0927
with timer("reading one probe"):
    os.chdir( folderpath1 )
    probefile1 = 'probe_00118.dat'
    Lsep1      = 13.12627
    probe1 = ProbeData( probefile1 )
#1125
with timer("reading one probe"):
    os.chdir( folderpath1 )
#    probefile1 = 'probe_00176.dat'
    probefile1 = 'probe_00189.dat'
    Lsep1      = 13.12627
    probe1 = ProbeData( probefile1 )
"""   
    
#1221
with timer("reading one probe"):
    os.chdir( folderpath1 )
#    probefile1 = 'probe_00176.dat'
    probefile1 = 'probe_00175.dat'
    Lsep1      = 13.266
    probe1 = ProbeData( probefile1 )    


with timer("clean data"):
    probe_s.cleandata(starttime=20)
    probe1.cleandata(starttime=20)
    
with timer('psd'):
    
    os.chdir(outpath1)
    
    probe_s.psd( 8, 0.5 )
    probe1.psd( 8, 0.5 )
    
    fig, ax = plt.subplots( figsize=[8,8], constrained_layout=True )
    
    St_Lsep_s = probe_s.St * Lsep_s
    St_Lsep1 = probe1.St * Lsep1
    
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

    ax.semilogx(St_Lsep_s, probe_s.nd_pprime_fwpsd,'red', linewidth=1)
    ax.semilogx(St_Lsep1, probe1.nd_pprime_fwpsd,'blue', linewidth=1)
    
    ax.set_xlim( [0.01,100] )
    ax.set_ylim( [0.0,1.0] )
    
    ax.set_xlabel( r'$St=f\cdot L_{sep}/U_{\infty}$', 
                  fontdict={'size':24} )
    ax.set_ylabel( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$', 
                  fontdict={'size':24} )
    if Pure:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        
    ax.grid(visible=True, which='both',axis='both',color='gray',
            linestyle='--',linewidth=0.2)
    
    plt.savefig("1221_psd")
    
    plt.show()
  
#    print(probe1.nperseg)
    
#    print(probe1.df['pprime'])

#    probe1.write_psd()

#%% fix probe data without headers of Luis case
'''
os.chdir( folderpath0 )

os.chdir( os.pardir )

filelist1 = get_filelist( folderpath1 )

headers = list()

for file in filelist1:
    
    with open( file, 'r' ) as f:
        
        headers.append( f.readline() )

filelist0 = get_filelist( folderpath0 )

# test with one file

#with open('probe_00097.dat','r+') as f:
#    old = f.read()
#    f.seek(0)
#    f.write(headers[0])
#    f.write(old)

for i, file in enumerate( filelist0 ):
    
    with open( file, 'r+' ) as f:
    
        old = f.read()
        f.seek(0)
        f.write(headers[i])
        f.write(old)
'''