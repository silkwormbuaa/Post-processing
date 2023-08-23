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

datapath = '/home/wencanwu/my_simulation/tutorials/tutorial_cylinder/results/'

gridfile = datapath + 'inca_grid.bin'

ibfile = datapath + 'cube_NEW.stl'

G = GridData( gridfile )

G.verbose = False

G.read_grid()

ibmsh = read_stl( ibfile )

x_init = [50.0, 50.0, 1.0]

bl_g = G.g[2]

print(bl_g.gx)
print(bl_g.gy)

bl_g.vol_fra = ray_tracer( ibmsh, bl_g, x_init, 1.0, verbose=True )

print( bl_g.vol_fra )

print(f"assign bl_g {bl_g.num}")