#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_yslice.py
@Time    :   2025/03/16 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   show contours of y-normal plane statistics
'''

import os
import sys
import numpy              as     np
import pandas             as     pd
import pyvista            as     pv
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.params      import Params
from   vista.directories import Directories
from   vista.statistic   import StatisticData

# =============================================================================
def main():
    
    casefolder  = '/home/wencan/temp/231124/'
    loc         = 0.1
    sub_folder  = 'y_' + str(loc)
    vars_in     = ['pp','p','u','v','w','T']
    bbox        = [-100,0,-2,10,-20,20]

    dirs       = Directories( casefolder )
    params     = Params( dirs.case_para_file )
    grid       = GridData( dirs.grid )
    
    grid.read_grid()
    blocklist, index = grid.select_sliced_blockgrids( 'Y', loc, bbox=bbox )
    
    os.chdir( dirs.pp_statistics + '/' + sub_folder )
    stat        = StatisticData( f'stat_{sub_folder}.bin' )
    stat.grid3d = grid
    stat.read_statistic( blocklist, vars_in=vars_in )
    
    dataset = pv.MultiBlock(stat.create_vtk_multiblock( blocklist, vars=vars_in ))
    
    visualize( dataset )
    

def visualize( dataset:pv.MultiBlock ):
    
    dataset = dataset.cell_data_to_point_data().combine()
    
    
    p = pv.Plotter()
    p.add_mesh( dataset, scalars='u', cmap='coolwarm' )

    p.add_axes()
    p.show()
    p.close()



# =============================================================================
if __name__ == "__main__":
    main()