#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   zslice_pkl.py
@Time    :   2024/02/28 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import pickle
import numpy             as     np
from   mpi4py            import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.tools       import get_filelist
from   vista.tools       import distribute_mpi_work
from   vista.plot_style  import plot_slicez_stat




        cbar = r'$u/u_{\infty}$'
        cbar_levels = np.linspace( -0.2, 1, 37)
        cbar_ticks  = np.linspace( -0.2, 1, 7)
        tag = f't = {snap.itime:6.2f} ms'
        plot_slicez_stat( xx,yy,u/507,
                          filename=f'u_{snap.itstep:08d}',
                          col_map='coolwarm',
                          cbar_label=cbar,
                          separation=sep_line_file,
                          sonic=False,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          tag=tag,
                          tag_loc=[-13,6],
                          x_lim=[-15,10],
                          y_lim=[0,8],
                          pure=False)
        
        cbar = r'$DS$'
        cbar_levels = np.linspace( 0.0, 0.8,33)
        
        plot_slicez_stat( xx,yy,DS,
                          filename=f'DS_{snap.itstep:08d}',
                          col_map='Greys_r',
                          cbar_label=cbar,
                          separation=sep_line_file,
                          sonic=False,
                          cbar_levels=cbar_levels,
                          tag=tag,
                          tag_loc=[-13,6],
                          x_lim=[-15,10],
                          y_lim=[0,8],
                          pure=False)
        
        
        os.system(f'mv u_{snap.itstep:08d}.png ./figures_u/')
        os.system(f'mv DS_{snap.itstep:08d}.png ./figures_DS/')

