#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   write_zslice.py
@Time    :   2024/10/28 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Read in a 3D statistics.bin file and write out mid-plane z-slice with
             all the statistics variables.
'''


import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.params      import Params
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder

# =============================================================================

def main():
    case_dirs = ['smooth_mid',      '231124','241030','241018',
                 'smooth_adiabatic','221014','220926','220825','220927','221221',
                 '240210',          '240211']
    vars_out  = ['u','v','w','p','T','mach','tke',
                 'u`u`','v`v`','w`w`','u`v`','p`']

    for case in case_dirs:
        write_incip_plane( '/home/wencan/temp/' + case, vars_out )
        print(f"Finished writing {case} incip plane data.\n")

# =============================================================================

def write_incip_plane( casedir, vars_out ):

    os.chdir( casedir )
    dirs    = Directories( casedir )
    params  = Params( dirs.case_para_file )
    x_imp   = params.x_imp
    delta   = params.delta_0
    rescale = [-x_imp, 0.0, 0.0, delta, delta, delta]
    x_incip = params.x_incip*delta + x_imp
    
    grid = GridData( dirs.grid )
    grid.read_grid()
    blocklist, index = grid.select_sliced_blockgrids( 'X', x_incip )

    stat3d = StatisticData( dirs.statistics )
    stat3d.grid3d = grid
    stat3d.read_statistic( blocklist, vars_in=stat3d.full_vars )

    stat3d.compute_vars( blocklist, vars_new=['mach','RS', 'p`'] )
    stat2d = stat3d.get_slice( 'X', x_incip )
    
    os.chdir( create_folder(dirs.stat_incip) )
    stat2d.write_vtm( 'stat_incip.vtm', vars=vars_out, block_list=blocklist,
                      rescale=rescale )

# =============================================================================
if __name__ == "__main__":
    main()