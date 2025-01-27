#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   new_grids.py
@Time    :   2025/01/27
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   create a new mesh for convergent-divergent riblet case 250127
'''


import os
import sys
import shutil
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.bloxx       import Grid_bloxx
from   vista.bloxx       import Mesh_bloxx
from   vista.directories import create_folder


# -- write a total new mesh

output_folder = '/home/wencanwu/my_simulation/STBLI_low_Re/250127/grid'

mesh = Mesh_bloxx()

domain_min = [-116, 0, -10.4]
domain_max = [ 120, 86, 10.4]

# build the bottom layer of blocks

lxs = np.linspace( domain_min[0], domain_max[0], 129 )
lzs = np.linspace( domain_min[2], domain_max[2], 9   )

for i in range(len(lxs)-1):
    for j in range(len(lzs)-1):
        
        grid = Grid_bloxx()
        grid.variables['LX'] = f"{lxs[i]:.8f}, {lxs[i+1]:.8f},"
        grid.variables['LY'] = '-1.25760000,   0.00000000,'
        grid.variables['LZ'] = f"{lzs[j]:.8f}, {lzs[j+1]:.8f},"

        mesh.grids.append(grid)

# build the second layer of blocks

lxs = np.linspace( domain_min[0], domain_max[0], 129 )
lzs = np.linspace( domain_min[2], domain_max[2], 9   )

for i in range(len(lxs)-1):
    for j in range(len(lzs)-1):
        
        grid = Grid_bloxx()
        grid.variables['LX'] = f"{lxs[i]:.8f}, {lxs[i+1]:.8f},"
        grid.variables['LY'] = '0.00000000, 1.7369534200,'
        grid.variables['LZ'] = f"{lzs[j]:.8f}, {lzs[j+1]:.8f},"

        mesh.grids.append(grid)

# build the third layer of blocks

lxs = np.linspace( domain_min[0], domain_max[0], 65 )
lzs = np.linspace( domain_min[2], domain_max[2], 9  )

for i in range(len(lxs)-1):
    for j in range(len(lzs)-1):
        
        grid = Grid_bloxx()
        grid.variables['LX'] = f"{lxs[i]:.8f}, {lxs[i+1]:.8f},"
        grid.variables['LY'] = '1.7369534200, 5.0103126500,'
        grid.variables['LZ'] = f"{lzs[j]:.8f}, {lzs[j+1]:.8f},"

        mesh.grids.append(grid)

# build the fourth layer of blocks

lxs = np.linspace( domain_min[0], domain_max[0], 33 )
lzs = np.linspace( domain_min[2], domain_max[2], 9  )

for i in range(len(lxs)-1):
    for j in range(len(lzs)-1):
        
        grid = Grid_bloxx()
        grid.variables['LX'] = f"{lxs[i]:.8f}, {lxs[i+1]:.8f},"
        grid.variables['LY'] = '5.0103126500, 11.1790910000,'
        grid.variables['LZ'] = f"{lzs[j]:.8f}, {lzs[j+1]:.8f},"
        
        mesh.grids.append(grid)

# build the fifth layer of blocks

lxs = np.linspace( domain_min[0], domain_max[0], 17 )
lzs = np.linspace( domain_min[2], domain_max[2], 9  )

for i in range(len(lxs)-1):
    for j in range(len(lzs)-1):
        
        grid = Grid_bloxx()
        grid.variables['LX'] = f"{lxs[i]:.8f}, {lxs[i+1]:.8f},"
        grid.variables['LY'] = '11.1790910000, 22.8044042000,'
        grid.variables['LZ'] = f"{lzs[j]:.8f}, {lzs[j+1]:.8f},"
        
        mesh.grids.append(grid)

# build the sixth layer of blocks

lxs = np.linspace( domain_min[0], domain_max[0], 17 )
lzs = np.linspace( domain_min[2], domain_max[2], 5  )

for i in range(len(lxs)-1):
    for j in range(len(lzs)-1):
        
        grid = Grid_bloxx()
        grid.variables['LX'] = f"{lxs[i]:.8f}, {lxs[i+1]:.8f},"
        grid.variables['LY'] = '22.8044042000, 44.7127788000'
        grid.variables['LZ'] = f"{lzs[j]:.8f}, {lzs[j+1]:.8f},"
        
        mesh.grids.append(grid)
        
# build the seventh layer of blocks

lxs = np.linspace( domain_min[0], domain_max[0], 17 )
lzs = np.linspace( domain_min[2], domain_max[2], 3  )

for i in range(len(lxs)-1):
    for j in range(len(lzs)-1):
        
        grid = Grid_bloxx()
        grid.variables['LX'] = f"{lxs[i]:.8f}, {lxs[i+1]:.8f},"
        grid.variables['LY'] = '44.7127788000, 86.0000000000'
        grid.variables['LZ'] = f"{lzs[j]:.8f}, {lzs[j+1]:.8f},"
        
        mesh.grids.append(grid)
        

# set the number of cells in each direction

for grid in mesh.grids:
    
    grid.variables['NX']     = '32,'
    grid.variables['NY']     = '32,'
    grid.variables['NZ']     = '32,'
    
    if grid.LY[0] >= 0.0:
        grid.variables['SHAPEY'] = '"LIN" ,'
        grid.variables['PARAMY'] = '1.02000000e+00  , 0.00000000e+00  , 0.00000000e+00  ,'

# set boundary conditions

    # - inlet

    if grid.LX[0] == domain_min[0]:
        
        if grid.LY[0] >= 10.0:
            grid.variables['BX1'] = '"RI_INFLOW" ,'
        elif grid.LY[1] <= 0.0:
            grid.variables['BX1'] = '"DUMMY" ,'
        else:
            grid.variables['BX1'] = '"DF_INFLOW" ,'
    
    # - sides
    
    if grid.LZ[0] == domain_min[2]:
        grid.variables['BZ1'] = '"CYC" ,'
    
    if grid.LZ[1] == domain_max[2]:
        grid.variables['BZ2'] = '"CYC" ,'
    
    # - upper
    
    if grid.LY[1] == domain_max[1]:
        grid.variables['BY2'] = '"RI_INFLOW" ,'
        
    # - outlet
    
    if grid.LX[1] == domain_max[0]:
        grid.variables['BX2'] = '"OUTFLOW" ,'

mesh.sort_grids()
shutil.rmtree( output_folder )
mesh.save_grid( create_folder(output_folder) )

