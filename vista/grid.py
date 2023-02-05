#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   grid.py
@Time    :   2023/02/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import sys

import os

import numpy             as np

import pandas            as pd

sys.path.append('..')

from utils.read_binary   import read_int_bin

from utils.read_binary   import read_flt_bin

from utils.read_binary   import read_log_bin

from utils.read_binary   import read_chr_bin

# ----------------------------------------------------------------------
# >>> Initialize a grid from binary file                        ( 1-0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

class GridData:
    
    def __init__(self, grid_dir ):
        
        # directory to the grid binary file
        self.dir = grid_dir
        
        # grid file size
        self.fsize = os.stat( grid_dir ).st_size
        
        # file pointer position
        self.pos = 0
        
        # total number of grids with 3 layers ghost cells 
        self.n_total = 0
        
        # number of grids without ghost cells
        self.n_real  = 0
        
        # BlockGrid List
        self.g = list()
        
        # verbose ?
        self.verbose = False

# ----------------------------------------------------------------------
# >>> Read Grid File Header                                     ( 1-1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_grid_header( self, file ):
        
        sfl = 8
        sin = 4
        slg = 4 
        # sch = 1 by default
        
        # read file format
        # Since IO Fortran will write 4 extra spaces before and after 
        # each 'write', which indicates the length of output and also 
        # final check of length of output. So spaces should be skipped
        # when reading binary data.
        
        self.pos = 0
        
        self.pos += 4        
        file.seek( self.pos )
        
        self.grid_file_format = read_int_bin( file.read(sin), sin )
        self.pos += sin
        
        self.pos += 4
        self.pos += 4
        
        file.seek( self.pos )
        
        self.grid_with_solver = read_log_bin( file.read(slg), slg )
        self.pos += slg
        
        self.pos += 4
        self.pos += 4
        
        file.seek( self.pos )
        
        self.bl_num = read_int_bin( file.read(sin), sin )
        self.pos += sin
        
        self.pos += 4
        
        self.header_size = self.pos
        

        if self.verbose:
            
            print(self.grid_file_format)
            print(self.grid_with_solver)
            print(self.bl_num)

# ----------------------------------------------------------------------
# >>> Read Grid File Body                                       ( 1-2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/03  - created
#
# Desc
#
# ----------------------------------------------------------------------
    def read_grid_body( self, file ):
        
        end_of_file = False
        
        # new position after reading header
        
        file.seek( self.header_size )
        
        self.pos = self.header_size
        
        # give every blockgrid an index
        i = 0
        
#        self.g.append( BlockGrid( file, 0,self.grid_with_solver))
        
        while not end_of_file:
            
            # read in block grid one by one
            # input: file, index, if grid_with_solver
                        
            self.g.append( BlockGrid(file,i,self.grid_with_solver) )
            
            self.pos = self.pos + self.g[-1].size
            
            i += 1
            
            if self.pos >= self.fsize: end_of_file = True

# ----------------------------------------------------------------------
# >>> Get Sorted Block Grids Groups                             ( 1-3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/05  - created
#
# Desc
#
# - based on lx0,ly0,lz0, sorting block grids
# - same x-y location blocks form a block group
#
# ----------------------------------------------------------------------
    def get_sorted_groups( self ):
        
        if len(self.g) == 0:
            raise ValueError('Please read grid blocks first!')
        
        lx0 = list()
        ly0 = list()
        lz0 = list()
        bl_index = list()
        
        # get lx0,ly0,lz0 and bl_index lists
        for i in range(self.bl_num):
            
            lx0.append(self.g[i].lx0)
            ly0.append(self.g[i].ly0)
            lz0.append(self.g[i].lz0)
            bl_index.append(self.g[i].num)
        
        df = pd.DataFrame(lx0,columns=['lx0'])
        df['ly0'] = ly0
        df['lz0'] = lz0
        df['bl_index'] = bl_index
        
        # sort based on lx0,ly0,lz0
        df_sorted = df.sort_values(by=['lx0','ly0','lz0'])
        
        # group pd raws and aggregate lz0 and bl_index into lists
        grouped  = df_sorted.groupby(by=['lx0','ly0']).agg(list)
        
        grouped = grouped.reset_index()
        
        self.g_group = np.array( grouped['bl_index'] )
        

# ----------------------------------------------------------------------
# >>> Class Block Grid                                          ( 2-0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/03  - created
#
# Desc
#
# - initialize an object of block grid and read grid info
#
# ----------------------------------------------------------------------

class BlockGrid:
    
    def __init__( self, file , num , grid_with_solver ):
        
        self.verbose = False
        self.grid_with_solver = grid_with_solver
        
        sin = 4
        sfl = 8
        slg = 4
        
        # index of BlockGrid, like first block's grid, second ...
        self.num = num
        if self.verbose: print( 'grid number is', self.num )
        
        # size of this BlockGrid (in bytes)
        self.size = 0
        
        # pass spaces before 'write'
        file.read(4)
        self.size += 4
        
        # read grid dimension
        self.npx = read_int_bin( file.read(sin), sin )
        self.npy = read_int_bin( file.read(sin), sin )
        self.npz = read_int_bin( file.read(sin), sin )
        
        self.size += 3*sin
        if self.verbose: print( 'npx = ', self.npx )
        
        # read grid bounding box corner location
        self.lx0 = read_flt_bin( file.read(sfl), sfl )
        self.lx1 = read_flt_bin( file.read(sfl), sfl )
        self.ly0 = read_flt_bin( file.read(sfl), sfl )
        self.ly1 = read_flt_bin( file.read(sfl), sfl )
        self.lz0 = read_flt_bin( file.read(sfl), sfl )
        self.lz1 = read_flt_bin( file.read(sfl), sfl )
        
        self.size += 6*sfl
        
        # read parameters of grid: 3 components in each direction
        self.paramx = read_flt_bin( file.read(3*sfl), sfl )
        self.paramy = read_flt_bin( file.read(3*sfl), sfl )
        self.paramz = read_flt_bin( file.read(3*sfl), sfl )
        
        self.size += 9*sfl
        if self.verbose: print( 'paramx = ', self.paramx )
        
        # read shapes of grids: char(len=10) for each direction
        self.shapex = read_chr_bin( file.read(10) )
        self.shapey = read_chr_bin( file.read(10) )
        self.shapez = read_chr_bin( file.read(10) )
        
        self.size += 30
        if self.verbose: print( 'shapex = ', self.shapex )
        
        # read boundary values: char(len=10)*2 for each direction
        self.bx1 = read_chr_bin( file.read(10) )
        self.bx2 = read_chr_bin( file.read(10) )
        self.by1 = read_chr_bin( file.read(10) )
        self.by2 = read_chr_bin( file.read(10) )
        self.bz1 = read_chr_bin( file.read(10) )
        self.bz2 = read_chr_bin( file.read(10) )
        
        self.size += 60
        if self.verbose: print( 'bx1 = ', self.bx1 )
        
        # read boundary parameters: int*3 for each one
        self.px1 = read_int_bin( file.read(3*sin), sin )
        self.px2 = read_int_bin( file.read(3*sin), sin )
        self.py1 = read_int_bin( file.read(3*sin), sin )
        self.py2 = read_int_bin( file.read(3*sin), sin )
        self.pz1 = read_int_bin( file.read(3*sin), sin )
        self.pz2 = read_int_bin( file.read(3*sin), sin )
        
        self.size += sin*3*6
        if self.verbose: print( 'px1 = ', self.px1 )
        
        # read times_pi (if lx * pi ?)
        self.times_pi = read_log_bin( file.read(slg), slg )
        
        self.size += slg
        if self.verbose: print( 'times_pi = ', self.times_pi )
        
        # read boundary values: float*7 for each one
        self.valx1 = read_flt_bin( file.read(7*sfl), sfl )
        self.valx2 = read_flt_bin( file.read(7*sfl), sfl )
        self.valy1 = read_flt_bin( file.read(7*sfl), sfl )
        self.valy2 = read_flt_bin( file.read(7*sfl), sfl )
        self.valz1 = read_flt_bin( file.read(7*sfl), sfl )
        self.valz2 = read_flt_bin( file.read(7*sfl), sfl )
        
        self.size += sfl*7*6
        if self.verbose: print( 'valx1 = ', self.valx1 )
        
        # read transient: logic for each one
        self.transx1 = read_log_bin( file.read(slg), slg )
        self.transx2 = read_log_bin( file.read(slg), slg )
        self.transy1 = read_log_bin( file.read(slg), slg )
        self.transy2 = read_log_bin( file.read(slg), slg )
        self.transz1 = read_log_bin( file.read(slg), slg )
        self.transz2 = read_log_bin( file.read(slg), slg )
        
        self.size += slg*6
        if self.verbose: print( 'transx1 = ', self.transx1 )
        
        # read fluid names: char(len=15) for each one
        self.fluidx1 = read_chr_bin( file.read(15) )
        self.fluidx2 = read_chr_bin( file.read(15) )
        self.fluidy1 = read_chr_bin( file.read(15) )
        self.fluidy2 = read_chr_bin( file.read(15) )
        self.fluidz1 = read_chr_bin( file.read(15) )
        self.fluidz2 = read_chr_bin( file.read(15) )
        
        self.size += 90
        if self.verbose: print( 'fluidx1 = ', self.fluidx1 )
        
        # read extra space after 'write'
        file.read(4)
        self.size += 4    
        
        
        # check if grid_with_solver, if so, extra info should be read
        # for current setup, no need to read.
        if self.grid_with_solver:
            pass
        
        
        # read extra space before 'write'
        file.read(4)
        self.size += 4
        
        # read grid geometric parameters
        # geometric parameters include ghost cells
        # - gx(y,z) cell center coordinates
        # - dx(y,z) distance between cell centers
        # - hx(y,z) cell widths
        # - wx(y,z) linear interpolation coefficients
        
        self.gx = read_flt_bin( file.read( (self.npx+6)*sfl ), sfl )
        self.dx = read_flt_bin( file.read( (self.npx+6)*sfl ), sfl )
        self.hx = read_flt_bin( file.read( (self.npx+6)*sfl ), sfl )
        self.wx = read_flt_bin( file.read( (self.npx+6)*sfl ), sfl )
        
        self.gy = read_flt_bin( file.read( (self.npy+6)*sfl ), sfl )
        self.dy = read_flt_bin( file.read( (self.npy+6)*sfl ), sfl )
        self.hy = read_flt_bin( file.read( (self.npy+6)*sfl ), sfl )
        self.wy = read_flt_bin( file.read( (self.npy+6)*sfl ), sfl )

        self.gz = read_flt_bin( file.read( (self.npz+6)*sfl ), sfl )
        self.dz = read_flt_bin( file.read( (self.npz+6)*sfl ), sfl )
        self.hz = read_flt_bin( file.read( (self.npz+6)*sfl ), sfl )
        self.wz = read_flt_bin( file.read( (self.npz+6)*sfl ), sfl )
        
        self.size += ( self.npx+self.npy+self.npz + 18 ) * sfl * 4
        if self.verbose: print( 'self.gy = ', self.gy )
        
        # read extra spaces after 'write'
        file.read(4)
        self.size += 4
        
        