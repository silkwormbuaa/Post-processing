#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   forces.py
@Time    :   2025/02/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   plot streamwise distribution of pressure and friction forces 
             (data from forces on the IB of each block)
'''

import os
import sys
import pickle
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.timer       import timer
from   vista.params      import Params
from   vista.grid        import GridData
from   vista.line        import LineData
from   vista.directories import Directories
from   vista.tools       import get_filelist
from   vista.directories import create_folder

class Bl_Force:
    
    def __init__(self, filename):
        
        with open(filename, 'r') as f:
            head_line = f.readline().strip()
            headers = None    
            if head_line[0] == '#':
                headers = head_line.lstrip('#').split()
        
        if not headers:
            headers = ['time', 'fric_x', 'fric_y', 'fric_z', 'pres_x', 'pres_y', 'pres_z']
            self.df = pd.read_csv( filename, delimiter=r'\s+', names=headers )
        
        else:
            self.df = pd.read_csv( filename, delimiter=r'\s+', names=headers, skiprows=1)
            
        self.ave_fric_x = self.df['fric_x'].mean()
        self.ave_fric_y = self.df['fric_y'].mean()
        self.ave_fric_z = self.df['fric_z'].mean()
        self.ave_pres_x = self.df['pres_x'].mean()
        self.ave_pres_y = self.df['pres_y'].mean()
        self.ave_pres_z = self.df['pres_z'].mean()
        
        self.index = int( os.path.basename(filename).split('_')[-1].split('.')[0] )
        
        print(f'read in {filename} {self.index}.')

# =============================================================================

casefolder    = '/home/wencan/temp/250120'
dirs          = Directories( casefolder )
forces_files  = get_filelist( dirs.for_dir, 'force_0' )
params        = Params( dirs.case_para_file )
dynp          = 0.5*params.rho_ref*params.u_ref**2

datas = '/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/wall_projection/streamwise_vars.pkl'
lines = LineData()
with open( datas, 'rb' ) as f:    lines.df = pickle.load( f )

os.chdir( create_folder(dirs.pp_forces) )

if not os.path.exists('mean_forces_x.csv'):

    grid      = GridData( dirs.grid )
    grid.read_grid()

    block_list, _ = grid.select_sliced_blockgrids( 'Y', 0.0 )

    grouped_block_list = grid.group_by_range('xy', block_list=block_list )

    colblock_list = list(list(zip(*grouped_block_list))[0])

    x = [0.5*(grid.g[num-1].lx0+grid.g[num-1].lx1) for num in colblock_list ]


    bl_forces = [ Bl_Force(file) for file in forces_files ]

    index_bl_forces = [bl_force.index for bl_force in bl_forces]

    columns = ['x', 'mean_fric_x', 'mean_fric_y', 'mean_fric_z', 'mean_pres_x', 'mean_pres_y', 'mean_pres_z']
    df = pd.DataFrame(columns=columns)

    for i, rowblock_list in enumerate(grouped_block_list):
        mean_fric_x = np.mean([bl_forces[index_bl_forces.index(num)].ave_fric_x for num in rowblock_list])
        mean_fric_y = np.mean([bl_forces[index_bl_forces.index(num)].ave_fric_y for num in rowblock_list])
        mean_fric_z = np.mean([bl_forces[index_bl_forces.index(num)].ave_fric_z for num in rowblock_list])
        mean_pres_x = np.mean([bl_forces[index_bl_forces.index(num)].ave_pres_x for num in rowblock_list])
        mean_pres_y = np.mean([bl_forces[index_bl_forces.index(num)].ave_pres_y for num in rowblock_list])
        mean_pres_z = np.mean([bl_forces[index_bl_forces.index(num)].ave_pres_z for num in rowblock_list])
        newrow = [x[i], mean_fric_x, mean_fric_y, mean_fric_z, mean_pres_x, mean_pres_y, mean_pres_z]
        newrow_df = pd.DataFrame([newrow], columns=columns)
        df = pd.concat( [df, newrow_df], ignore_index=True )

    g = grid.g[colblock_list[0]-1]
    block_area = (g.lx1 - g.lx0) * (g.lz1 - g.lz0)

    for var in df.columns[1:]:
        df[var] = df[var] / block_area

    df.to_csv( 'mean_forces_x.csv', sep='\t',index=False )
    print(f'Save mean_forces_x.csv in {os.getcwd()}.')

else:
    df = pd.read_csv('mean_forces_x.csv', delimiter=r'\s+')
    print(f'Read mean_forces_x.csv in {os.getcwd()}.')
    
print( df )

# ----------------------------------------------------------------------
# >>> total force in x                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/02/13  - created
#
# Desc
#
# ----------------------------------------------------------------------
plt.rcParams['font.size']           = 20

total_force_x = df['mean_fric_x'] + df['mean_pres_x']

fig, ax = plt.subplots()

x = (df['x']-50.4)/5.2

ax.plot( x, df['mean_fric_x']/dynp*1000, label='friction force in x', color='red' )
ax.plot( x, df['mean_pres_x']/dynp*1000, label='pressure force in x', color='blue' )
ax.plot( x, total_force_x/dynp*1000,     label='total force in x', color='black' ,marker='o')
ax.hlines(0.0, x.min(), x.max(), color='gray', linestyle='--')

# - smooth wall

ax.plot( lines.df['x'], lines.df['Cf'], label='smooth wall fric', color='gray', linestyle='-' )


ax.set_xlabel(r'$(x-x_\delta)/\delta_0$')
ax.set_ylabel(r'$\tau$')

# ax.set_xlim([-20, 10])
ax.set_ylim([-20,40])


plt.legend()
plt.show()
