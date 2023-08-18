#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   read_statistics.py
@Time    :   2023/08/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid               import GridData

from   vista.geometry           import read_stl

from   vista.geometry           import ray_tracer

datapath = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/'

gridfile = datapath + 'inca_grid.bin'

ibfile = datapath + 'ib.stl'

G = GridData( gridfile )

G.verbose = True

G.read_grid()

ibmsh = read_stl( ibfile )

x_init = [0.0, 10.0, 0.0]

bl_g = G.g[0]

bl_g.vol_fra = ray_tracer( ibmsh, bl_g, x_init, 1.0,verbose=True )

print(f"assign bl_g {bl_g.num}")