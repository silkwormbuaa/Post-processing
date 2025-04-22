#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   check_skewness.py
@Time    :   2025/04/22 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   want to check the skewness of the velocity signals from probes
             expect to see near the shock foot, the skewness changes side across
             the shock foot.
'''


import os
import sys
import numpy             as     np
from   scipy.stats       import skew, kurtosis

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.probe       import ProbeData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.directories import create_folder


def main():

    # =============================================================================

    casefolder = '/home/wencan/temp/241030'
    vars       = ['u','v','w','rho', 'p', 'T']

    # =============================================================================

    dirs = Directories( casefolder )

    params = Params( dirs.case_para_file )
    withT  = params.prb_withT
    
    pfmax_i = params.prb_pfmax

    # -- get all probe files
    prb_files = get_filelist( dirs.prb_dir, 'probe_' )
    prb_files = prb_files[ pfmax_i-10:pfmax_i+10 ]
    
    os.chdir( create_folder( dirs.pp_signals ) )

    # ----- first get the index of probes signals at snapshot time points.

    prb_data = ProbeData( prb_files[0], withT=withT )
    prb_data.cleandata( t_start=20.0 )
    index    = prb_data.time_index( np.linspace(20.0, 61.0, 4101) )

    clock    = timer("Plot signals from probes")
    for k, file in enumerate(prb_files):
        
        prb    = ProbeData( file, withT=withT )
        prb.cleandata( 20.0 )
        prb.df = prb.df.iloc[index]
        
        skewness = skew( prb.df['p'] )
        print(f"skewness of {file} is {skewness:.3f}")


        progress = (k+1)/len(prb_files)
        print(f"{k+1}/{len(prb_files) } is done."+ clock.remainder(progress))
        print("---------------------\n")
        sys.stdout.flush()


# =============================================================================
if __name__ == "__main__":

    main()

