# -*- coding: utf-8 -*-
'''
@File    :   block.py
@Time    :   2023/02/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Class block
'''

import numpy             as     np

import pandas            as     pd

from   .io_binary        import read_int_bin
from   .io_binary        import read_flt_bin
from   .io_binary        import read_3Dflt_bin
from   .io_binary        import read_log_bin


# ----------------------------------------------------------------------
# >>> Class Statistic Block Data                                ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/01  - created
#
# Desc
#
# - initialize an object of block data.
# - for selected block, block data is read and RS is calculated.
#
# ----------------------------------------------------------------------

class BlockData:
    
    # if verbose
    
    verbose = False
    
    # size of type 
    
    sin = 4
    sfl = 8
    slg = 4
    
    def __init__( self, file, block_list, n_var, vars, var_indx ):
        
        """
        file       : opened file object
        block_list : list of blocks's numbers
        n_var      : number vars in the stored blocks
        vars       : list of selected variable name strings
        var_indx   : list of indexes of selected variables in data chunk
        """
        
        # start position of file pointer
        
        pos_start = file.tell()
        
        # index of block, can be read from blockdata itself
        
        self.num = 0
        
        # size of this BlockData (in bytes)
        
        self.size = 0
        self.n_var = n_var
        
        # empty list for future use
        
        self.df = pd.DataFrame(columns=vars)
        
        # matrix of friction projection on x-z plane

        self.df_fric = None
        
        # read global block number and block dimensions
        
        self.num = read_int_bin( file.read(self.sin), self.sin )
        self.npx = read_int_bin( file.read(self.sin), self.sin )
        self.npy = read_int_bin( file.read(self.sin), self.sin )
        self.npz = read_int_bin( file.read(self.sin), self.sin )
        
        self.np  = self.npx * self.npy * self.npz
        
        # if this block chunk will be read?
        # by default, to_fill is False
        self.to_fill = False
        if self.num in block_list: self.to_fill = True

        # primitive variables 'u,v,w,rho,rhoE'
        # mean variables = primitive variables + 'p T mu'
        if self.to_fill:
            
            for i, var in enumerate(vars):
                
                pos = pos_start + 4*self.sin + var_indx[i]*self.np*self.sfl
                
                file.seek( pos )
                
                buff = read_flt_bin( file.read(self.np*self.sfl), self.sfl )

                self.df[var] = buff
                
            
        # calculate the block data size in byte
        
        self.size = self.np*self.n_var*self.sfl + 4*self.sin
        
        # move file pointer to the end of current block

        file.seek( pos_start + self.size)


# ----------------------------------------------------------------------
# >>> Define a subclass SnapBlock                                 (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/12  - created
#
# Desc
#
# ----------------------------------------------------------------------

class SnapBlock(BlockData):
    
    def __init__( self, file, block_list, n_vars, vars, snap_with_gx, type, 
                  kind=4 ):

        """
        file         : opened file object
        block_list   : list of blocks's numbers
        n_vars       : number vars in the stored blocks
        vars         : list of selected variable name strings
        snap_with_gx : if contain grids
        type         : snapshot type, 'block' or 'slice'
        kind         : should be consistent to Snapshot.kind
        
        ! SnapBlock read in all variables.
        """
        
# ----- Initialize attributes of instance.

        # floating point precision format
        
        self.kind = kind

        # start position of file pointer
        
        pos_start = file.tell()
        
        self.pos_start = pos_start
        
        # start position of variable start
        
        self.pos_var_start = None
        
        # index of block, can be read from blockdata itself
        
        self.num = 0
        
        # size of this BlockData (in bytes)
        
        self.size = 0
        self.n_vars = n_vars
        
        # empty grids points list
        
        self.gx = None
        self.gy = None
        self.gz = None
        
        # empty dataframe for datachunk

        self.df = pd.DataFrame(columns=vars)
        
        # matrix of friction projection on x-z plane

        self.df_fric = None
        
# ----- read global block number and block dimensions
        
        self.num = read_int_bin( file.read(self.sin), self.sin )
        
        if self.num == 0:
            raise ValueError(f"Found a block with num==0 at {pos_start}.")
        
        dim = read_int_bin( file.read(3*self.sin), self.sin )
        self.npx = dim[0]
        self.npy = dim[1]
        self.npz = dim[2]
        
        self.np  = self.npx * self.npy * self.npz

        pos = pos_start + 4*kind
        
        if self.verbose:
            print(f"read in block {self.num}, dimension {dim}.")
        
# ----- read grid points

        n_grid = 0

        if snap_with_gx and type == 'block':
            
            self.gx = read_3Dflt_bin( pos, file, self.npx, 1, 1, kind )
            pos += self.npx*kind
            
            self.gy = read_3Dflt_bin( pos, file, 1, self.npy, 1, kind )
            pos += self.npy*kind
            
            self.gz = read_3Dflt_bin( pos, file, 1, 1, self.npz, kind )
            pos += self.npz*kind
            
            n_grid = self.npx + self.npy + self.npz
        
        if snap_with_gx and type == 'slice':
            
            # slic_type == 'X'
            
            if self.npx == 1: 
                
                self.gy = read_3Dflt_bin( pos, file, 1,self.npy,1, kind )
                pos += self.npy*kind

                self.gz = read_3Dflt_bin( pos, file, 1,1,self.npz, kind )
                pos += self.npz*kind
                
                n_grid = self.npy + self.npz

            # slic_type == 'Y' or 'W'
            
            elif self.npy == 1:

                self.gx = read_3Dflt_bin( pos, file, self.npx,1,1, kind )
                pos += self.npx*kind

                self.gz = read_3Dflt_bin( pos, file, 1,1,self.npz, kind )
                pos += self.npz*kind
                
                n_grid = self.npx + self.npz
            
            # slic_type == 'Z'
            
            elif self.npz == 1:
                
                self.gx = read_3Dflt_bin( pos, file, self.npx,1,1,kind )
                pos += self.npx*kind

                self.gy = read_3Dflt_bin( pos, file, 1,self.npy,1,kind )
                pos += self.npy*kind
                
                n_grid = self.npx + self.npy

# ----- record the starting position of variable chunk
            
        self.pos_var_start = pos
                
# ----- read variable data chunk

        # if this block chunk will be read?
        # by default, to_fill is False
        self.to_fill = False
        
        if (block_list is None) or (self.num in block_list):
            self.to_fill = True

        # primitive variables 'u,v,w,rho,rhoE'
        # mean variables = primitive variables + 'p T mu'
        if self.to_fill:
            
            for var in enumerate(vars):
            
                buff = read_flt_bin( file.read(self.np*kind), kind )

                self.df[var] = buff
                
# ----- recording and move file pointer
    
        # calculate the block data size in byte
        
        self.size = 4*self.sin + n_grid*kind + self.np*n_vars*kind
        
        # move file pointer to the end of current block(if not fill)

        file.seek( pos_start + self.size )
        