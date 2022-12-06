#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run4.py
@Time    :   2022/11/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Scripts for getting PSD
'''

import sys

sys.path.append('..')

import os

import matplotlib.pyplot as     plt

from   vista.psd         import ProbeData

from   utils.timer       import timer

from   utils.tools       import get_filelist


folderpath0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/probes/probe_x'
folderpath1 = '/home/wencanwu/my_simulation/temp/221125_lowRe/probes/probe_x'

outpath0 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/probes/psd_x_new'
outpath1 = '/home/wencanwu/my_simulation/temp/221125_lowRe/probes/psd_x'
#%%

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
            
            probe.psd( 4, 0.5 )

            probe.write_psd()

#%% print single line
"""
#
with timer("reading smooth wall probe"):
    os.chdir( folderpath0 )
    probefile_s  = 'probe_00142.dat'
    Lsep_s       = 13.12627403
    probe_s = ProbeData( probefile_s )

#1014
#with timer("reading one probe"):
#    os.chdir( folderpath1 )
#    probefile1 = 'probe_00144.dat'
#    Lsep1      = 9.805522
#    probe1 = ProbeData( probefile1 )

#0926
#with timer("reading one probe"):
#    os.chdir( folderpath1 )
#    probefile1 = 'probe_00145.dat'
#    Lsep1      = 9.676337
#    probe1 = ProbeData( probefile1 )

#0825
#with timer("reading one probe"):
#    os.chdir( folderpath1 )
#    probefile1 = 'probe_00135.dat'
#    Lsep1      = 11.340638
#    probe1 = ProbeData( probefile1 )

#0927
#with timer("reading one probe"):
#    os.chdir( folderpath1 )
#    probefile1 = 'probe_00118.dat'
#    Lsep1      = 13.12627
#    probe1 = ProbeData( probefile1 )
    
#1125
with timer("reading one probe"):
    os.chdir( folderpath1 )
    probefile1 = 'probe_00118.dat'
    Lsep1      = 13.12627
    probe1 = ProbeData( probefile1 )
    
with timer("clean data"):
    probe_s.cleandata(starttime=20)
    probe1.cleandata(starttime=20)
    
with timer('psd'):
    
    os.chdir(os.pardir)
    
    probe_s.psd( 8, 0.5 )
    probe1.psd( 8, 0.5 )
    
    fig, ax = plt.subplots( figsize=[10,8], constrained_layout=True )
    
    St_Lsep_s = probe_s.St * Lsep_s
    St_Lsep1 = probe1.St * Lsep1
    
    ax.minorticks_on()

    ax.semilogx(St_Lsep_s, probe_s.nd_pprime_fwpsd,'r', linewidth=1)
    ax.semilogx(St_Lsep1, probe1.nd_pprime_fwpsd,'b', linewidth=1)
    
    ax.set_xlim( [0.001,100] )
    
    ax.set_xlabel( r'$St=f\cdot L_{sep}/U_{\infty}$', 
                  fontdict={'size':24} )
    ax.set_ylabel( r'$f \cdot PSD(f)/ \int PSD(f) \mathrm{d} f$', 
                  fontdict={'size':24} )
    
    ax.grid( True, which = 'both', ls='--' )
    
    plt.savefig("0927_psd_xsep")
    
    plt.show()
  
#    print(probe1.nperseg)
    
#    print(probe1.df['pprime'])

#    probe1.write_psd()
"""
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