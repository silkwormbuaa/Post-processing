#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   cross_corelation_analysis.py
@Time    :   2024/09/09 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Cross correlation analysis of the shock motion and bubble size
'''

import os
import sys
import pickle
import numpy              as     np
import pandas             as     pd
import matplotlib.pyplot  as     plt
from   scipy.signal       import correlate

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.params       import Params
from   vista.probe        import ProbeData
from   vista.directories  import Directories
from   vista.math_opr     import find_parabola_max
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(fontsize=30)

def main():
    case_dir = '/home/wencan/temp/241030/'
    dirs     = Directories( case_dir )
    params   = Params( dirs.case_para_file )

    df       = read_prb( dirs.fetch_prb_from_type('pfmax'), params )
    plot_flucs( df['time'], [df['p_fluc'],df['u_fluc']], ['p','u'], 'pu_pfmax' )


def plot_flucs( times, flucs, vars, title='' ):
    
    fig, ax = plt.subplots( figsize=(15,8) )

    for i, fluc in enumerate(flucs):
        norm = np.std( fluc )
        ax.plot( times, np.array(fluc)/norm, label=vars[i] )

    ax.set_ylim( -2.1, 2.1 )
    ax.set_title(title)
    ax.legend()

    plt.show()
    plt.close()
    
    

def read_prb( file, params: Params ):
    
    # - read in the probe data 
    
    prb_data = ProbeData( file, params.prb_withT )
    prb_data.cleandata( t_start=20.0 )
    prb_data.get_fluc( prb_data.df.columns )
    index = prb_data.time_index( np.linspace(20.0, 61.0, 4101) )
    df    = prb_data.df.iloc[index]

    return df



# =============================================================================
if __name__ == "__main__":
    main()