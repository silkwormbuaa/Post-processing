#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   prb_signal_shockfoot.py
@Time    :   2025/08/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.params      import Params
from   vista.directories import Directories

def main():
    casefolder   = '/home/wencan/temp/250710'
    div_prb_file = 'probe_00158.dat'
    con_prb_file = 'probe_00583.dat'
    
    dirs = Directories( casefolder )
    params = Params( dirs.case_para_file )
    withT  = params.prb_withT
    
    div_prb = ProbeData( os.path.join(dirs.prb_dir, div_prb_file), withT=withT )
    con_prb = ProbeData( os.path.join(dirs.prb_dir, con_prb_file), withT=withT )

    div_prb.cleandata( t_start=20.0 )
    con_prb.cleandata( t_start=20.0 )
    
    # -- plot velocity signals
    
    fig, ax = plt.subplots(figsize=(12,8), constrained_layout=True)
    
    ax.plot( div_prb.df['time'], div_prb.df['u']/params.u_ref, 'black' )
    ax.plot( con_prb.df['time'], con_prb.df['u']/params.u_ref, 'red' )
    
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel(r"$u/u_{ref}$")
    
    plt.show()    
    plt.close()
    
    # -- plot pressure signals
    fig, ax = plt.subplots(figsize=(12,8), constrained_layout=True)
    ax.plot( div_prb.df['time'], div_prb.df['p']/params.p_ref, 'black' )
    ax.plot( con_prb.df['time'], con_prb.df['p']/params.p_ref, 'red' )
    
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel(r"$p/p_{ref}$")

    plt.show()
    plt.close()
    
    # -- plot vertical velocity signals
    fig, ax = plt.subplots(figsize=(12,8), constrained_layout=True)
    ax.plot( div_prb.df['time'], div_prb.df['v']/params.u_ref, 'black' )
    ax.plot( con_prb.df['time'], con_prb.df['v']/params.u_ref, 'red' )
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel(r"$v/u_{ref}$")
    plt.show()
    plt.close()
    
    # -- plot spanwisei velocity signals
    fig, ax = plt.subplots(figsize=(12,8), constrained_layout=True)
    ax.plot( div_prb.df['time'], div_prb.df['w']/params.u_ref, 'black' )
    ax.plot( con_prb.df['time'], con_prb.df['w']/params.u_ref, 'red' )
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel(r"$w/u_{ref}$")
    plt.show()
    plt.close()
    
def analyse_signal( df:pd.DataFrame ):
    
    '''
    Analyse the signal from a probe data frame.
    '''
    u = np.array( df['u'] )
    v = np.array( df['v'] )
    w = np.array( df['w'] )
    p = np.array( df['p'] )
    
        
    



# =============================================================================
if __name__ == "__main__":

    main()

