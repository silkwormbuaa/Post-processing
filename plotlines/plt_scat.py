#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_scat.py
@Time    :   2022/11/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Plotting scatter with solidity
'''

import os

import numpy             as     np

import matplotlib        as     mpl

import matplotlib.pyplot as     plt

from   plt_tools         import *

plt.rcParams['font.family'] = "serif"

#plt.rcParams.update({'font.size': 18})

OutPath = "/home/wencanwu/my_simulation/temp/221014_lowRe/DataPost"

DataFile = "solidity"

plt_DU   =  False

plt_vbar =  True 


os.chdir(OutPath)

data = PlotDataframe(DataFile)


if plt_DU :
    
    fig, ax = plt.subplots(figsize=[10,8])


    
    ax.scatter( data.df['Lambda'], 
                data.df['DU'], 
                marker='s' )

    
    ax.minorticks_on()
    
    ax.set_xscale( "log" )
    ax.set_xlabel( r"$\Lambda$", fontdict={'size':24} )
    
    ax.set_xlim( [0.01,1.0] )
        
    ax.set_ylabel( r'$\Delta U^+$', fontdict={'size':24} )
    
    ax.grid(True, which = 'both', ls='--')
    ax.set_axisbelow(True)
    
    plt.savefig("DU.png")
    plt.show()
    
if plt_vbar :
    
    fig, ax = plt.subplots(figsize=[10,8], constrained_layout=True )


    
    ax.scatter( data.df['D/δ'], 
                data.df['v_']*100, 
                marker='s' ,
                s = 100)

    
    ax.minorticks_on()
    
#    ax.set_xscale( "log" )
    ax.set_xlabel( r"$D/\delta_0$", fontdict={'size':36} )
    
    ax.set_xlim( [0.0,2.5] )
    ax.tick_params(axis='x',labelsize=32)
        
    ax.set_ylabel( r"${v}_{max}/u_{\infty}\times 10^2$", fontdict={'size':36} )
    ax.tick_params(axis='y',labelsize=32)
    
    
    ax.grid(True, which = 'major', ls='--')
    ax.set_axisbelow(True)
    
    plt.savefig("v_.png")
    plt.show()