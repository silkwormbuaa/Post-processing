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

from   .grid             import GridBlock


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
    
    # size of type 
    
    sin = 4
    sfl = 8
    slg = 4
    
    def __init__( self, file=None, block_list=None, n_var=None, vars=None, var_indx=None, verbose=False ):
        
        """
        can initialized as void or ...
        
        file       : opened file object
        block_list : list of blocks's numbers
        n_var      : number vars in the stored blocks
        vars       : list of selected variable name strings
        var_indx   : list of indexes of selected variables in data chunk
        """
        if file is None:
            pass
        else:
            self.init_from_file( file, block_list, n_var, vars, var_indx, verbose )
    
    
    def init_from_file( self, file, block_list, n_var, vars, var_indx, verbose ):
        
        # if verbose?
        self.verbose = verbose
        
        # start position of file pointer
        
        pos_start = file.tell()
        
        # size of this BlockData (in bytes)
        
        self.size = 0
        
        # number of variables
        
        self.n_var = n_var
        
        # empty list for future use
        
        self.df = pd.DataFrame(columns=vars)
        
        # matrix of friction projection on x-z plane

        self.df_fric = None
        
        # read global block number(index) and block dimensions
        
        self.num = read_int_bin( file.read(self.sin), self.sin )
        self.npx = read_int_bin( file.read(self.sin), self.sin )
        self.npy = read_int_bin( file.read(self.sin), self.sin )
        self.npz = read_int_bin( file.read(self.sin), self.sin )
        
        self.np  = self.npx * self.npy * self.npz
        
        # if this block chunk will be read?
        # by default, to_fill is False
        self.to_fill = False
        if self.num in block_list: self.to_fill = True
        
        if self.verbose:
            print(f"read in block {self.num},",end='')
            print(f"dimension {self.npx} {self.npy} {self.npz}. {self.to_fill}")

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
        
        # init grid points
        
        self.g = GridBlock()


# ----------------------------------------------------------------------
# >>> Drop ghost cells                                       (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/02/22  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def drop_ghost( self, buff=3 ):
        
# ----- check if the block is filled

        if self.df.empty:
            raise ValueError(f"BlockData {self.num} is empty.")
        
# ----- check if the block is 3D or 2D, if 2D, which type

        if any(x==1 for x in (self.npx,self.npy,self.npz)):
            if self.npx == 1: block_type = 'X'
            if self.npy == 1: block_type = 'Y'
            if self.npz == 1: block_type = 'Z'
        else:
            block_type = 'block'

# ----- drop ghost cells

        # N1, N2, N3 are dimensions with ghost cells
        
        N1 = self.npx
        N2 = self.npy
        N3 = self.npz
        n_var = len( self.df.columns )
        
        # drop ghost cell coordinates
        
        if self.g is not None:
            
            if block_type == 'block':
                self.g.gx = self.g.gx[buff:-buff]
                self.g.gy = self.g.gy[buff:-buff]
                self.g.gz = self.g.gz[buff:-buff]
            elif block_type == 'X':
                self.g.gy = self.g.gy[buff:-buff]
                self.g.gz = self.g.gz[buff:-buff]
            elif block_type == 'Y':
                self.g.gx = self.g.gx[buff:-buff]
                self.g.gz = self.g.gz[buff:-buff]
            elif block_type == 'Z':
                self.g.gx = self.g.gx[buff:-buff]
                self.g.gy = self.g.gy[buff:-buff]

        # drop ghost cells
        # Notice the binary data storage order [n_var, Z, Y, X]
        # Reshape to chunk -> slice -> reshape to vector
        
        if block_type == 'block':
            
            Nx = N1 - buff*2
            Ny = N2 - buff*2
            Nz = N3 - buff*2
            
            sol_buff = (self.df.values).T.reshape( n_var, N3, N2, N1 )
            sol_buff = sol_buff[ :, buff:-buff, buff:-buff, buff:-buff ]

        elif block_type == 'X':
            
            Nx = 1
            Ny = N2 - buff*2
            Nz = N3 - buff*2
            
            sol_buff = (self.df.values).T.reshape( n_var, N3, N2 )
            sol_buff = sol_buff[ :, buff:-buff, buff:-buff ]

        elif block_type == 'Y':
            
            Nx = N1 - buff*2
            Ny = 1
            Nz = N3 - buff*2

            sol_buff = (self.df.values).T.reshape( n_var, N3, N1 )
            sol_buff = sol_buff[ :, buff:-buff, buff:-buff ]
            
        elif block_type == 'Z':
            Nx = N1 - buff*2
            Ny = N2 - buff*2
            Nz = 1
        
            sol_buff = (self.df.values).T.reshape( n_var, N2, N1 )
            sol_buff = sol_buff[ :, buff:-buff, buff:-buff ]
        
        # Reshape the matrix ( no matter 2D or 3D) to long vectors
                
        sol_buff = sol_buff.reshape(( n_var, Nx*Ny*Nz ))
        
        df = pd.DataFrame(sol_buff.T, columns=self.df.columns)
        
        # Append to snap_data_clean
        
        self.npx = Nx
        self.npy = Ny
        self.npz = Nz
        self.df = df
        

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
    
    verbose = False
    
    def __init__( self, file=None, block_list=None, n_var=None, vars=None, 
                        snap_with_gx = None, type=None, kind=4 ):
        
        """
        file         : opened file object
        block_list   : list of blocks's numbers
        n_var       : number vars in the stored blocks
        vars         : list of selected variable name strings
        snap_with_gx : if contain grids
        type         : snapshot type, 'block' or 'slice'
        kind         : should be consistent to Snapshot.kind
        
        ! if nothing is given, build a void SnapBlock(), then fill_with_data()
        ! SnapBlock read in all variables.
        """
        
        if file is None:
            self.g = GridBlock()
        
        else:
            self.init_from_file( file, block_list, n_var, vars, snap_with_gx,
                                 type, kind)
        

    def init_from_file( self, file, block_list, n_var, vars, snap_with_gx, 
                        type, kind=4):
        
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
        self.n_var = n_var
        
        # empty grids points list
        
        n_grid = 0
        self.g = GridBlock()
        
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

        if snap_with_gx and type == 'block':
            
            self.g.gx = read_3Dflt_bin( pos, file, self.npx, 1, 1, kind )
            pos += self.npx*kind
            
            self.g.gy = read_3Dflt_bin( pos, file, 1, self.npy, 1, kind )
            pos += self.npy*kind
            
            self.g.gz = read_3Dflt_bin( pos, file, 1, 1, self.npz, kind )
            pos += self.npz*kind
            
            n_grid = self.npx + self.npy + self.npz
        
        if snap_with_gx and type == 'slice':
            
            # slic_type == 'X'
            
            if self.npx == 1: 
                
                self.g.gx = None
                
                self.g.gy = read_3Dflt_bin( pos, file, 1,self.npy,1, kind )
                pos += self.npy*kind

                self.g.gz = read_3Dflt_bin( pos, file, 1,1,self.npz, kind )
                pos += self.npz*kind
                
                n_grid = self.npy + self.npz

            # slic_type == 'Y' or 'W'
            
            elif self.npy == 1:

                self.g.gx = read_3Dflt_bin( pos, file, self.npx,1,1, kind )
                pos += self.npx*kind
                
                self.g.gy = None

                self.g.gz = read_3Dflt_bin( pos, file, 1,1,self.npz, kind )
                pos += self.npz*kind
                
                n_grid = self.npx + self.npz
            
            # slic_type == 'Z'
            
            elif self.npz == 1:
                
                self.g.gx = read_3Dflt_bin( pos, file, self.npx,1,1,kind )
                pos += self.npx*kind

                self.g.gy = read_3Dflt_bin( pos, file, 1,self.npy,1,kind )
                pos += self.npy*kind
                
                self.g.gz = None
                
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
            
            for var in vars:
            
                buff = read_flt_bin( file.read(self.np*kind), kind )

                self.df[var] = buff
                
# ----- recording and move file pointer
    
        # calculate the block data size in byte
        
        self.size = 4*self.sin + n_grid*kind + self.np*n_var*kind
        
        # move file pointer to the end of current block(if not fill)

        file.seek( pos_start + self.size )


# ----------------------------------------------------------------------
# >>> fill_with_data                                              (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/16  - created
#
# Desc
#
#   - initialize SnapBlock by filling data directly
# ----------------------------------------------------------------------

    def fill_with_data( self, num, dims, GX, df, type, kind=4 ):
        
        """
        num:  block number
        dims: [npx,npy,npz]
        GX:   [gx,gy,gz] or 2D format
        df:   dataframe
        type: snapshot type, 'block' or 'slice'
        """
        
        # matrix of friction projection on x-z plane

        self.df_fric = None
        
        # remember the kind
        
        self.kind = kind
        
# ----- fill in data
        
        self.num = num
        
        self.npx = dims[0]
        self.npy = dims[1]
        self.npz = dims[2]
        
        # grid coordinates vectors
        
        if type == 'block':
            self.g.gx = GX[0]; self.g.gy = GX[1]; self.g.gz = GX[2]
        
        if type == 'slice':
            if   self.npx == 1: 
                self.g.gx = None;  self.g.gy = GX[0]; self.g.gz = GX[1]
            elif self.npy == 1: 
                self.g.gx = GX[0]; self.g.gy = None;  self.g.gz = GX[1]
            elif self.npz == 1: 
                self.g.gx = GX[0]; self.g.gy = GX[1]; self.g.gz = None
        
        # dataframe
        
        self.df = df
            