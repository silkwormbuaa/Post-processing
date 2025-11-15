#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   check_convergence.py
@Time    :   2025/10/29 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Checking if data converge by plotting the mean changes with time.
'''

import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.probe       import ProbeFile
from   vista.directories import Directories
from   vista.params      import Params

def main():

    case_folder = '/home/wencan/temp/241030_z2'
    probe_file  = 'probe_00118.dat'
   
    dirs        = Directories( case_folder )
    params      = Params( dirs.case_para_file )

    os.chdir( dirs.prb_dir )
    
# -- check if probe header is available
    
    probe       = ProbeData( probe_file, withT=params.prb_withT )
    probes      = ProbeFile( dirs.set_prb )
    index       = int(probe_file.split('_')[-1].split('.')[0])
    loc         = probes.probes[ index-1 ].xyz
#    probe.cleandata( t_start=66.0 )
    
    p     = np.array( probe.df['p'] )
    itime = np.array( probe.df['time'] )
    
    pmean =  np.zeros_like( p )
    pmean[0] = p[0]
    
    for i in range(1,len(p)):
        pmean[i] = pmean[i-1]*i/(i+1) + p[i]/(i+1)    
        
    print(f"{probe_file} is located at {loc}, the mean value changes with time as the figure shows:")
    
    fig, ax = plt.subplots()
    
    ax.plot(itime,pmean/params.p_ref)
    ax.set_title( f"Mean pressure change with time at {loc}" )
    
    plt.show()
    plt.close()
    
    
    


# =============================================================================
if __name__ == "__main__":

    main()

