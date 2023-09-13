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
    
    def __init__( self, file, n_var, block_list, vars, var_indx ):
        
        self.verbose = False
        
        # start position of file pointer
        
        pos_start = file.tell()
        
        # index of block, can be read from blockdata itself
        self.num = 0
        
        # size of this BlockData (in bytes)
        self.size = 0
        self.n_var = n_var
        self.dim = []
        
        # empty list for future use
        self.df = pd.DataFrame(columns=vars)

        # size of type 
        
        sin = 4
        sfl = 8
        slg = 4        
        
        # read global block number and block dimensions
        
        self.num = read_int_bin( file.read(sin), sin )
        self.dim = read_int_bin( file.read(3*sin), sin )
        
        npx = self.dim[0]
        npy = self.dim[1]
        npz = self.dim[2]
        
        self.np  = npx * npy * npz
        
        # if this block chunk will be read?
        # by default, to_fill is False
        self.to_fill = False
        if self.num in block_list: self.to_fill = True

        # primitive variables 'u,v,w,rho,rhoE'
        # read primitive variables + 'p T mu', when fill is true.
        if self.to_fill:
            
            for i, var in enumerate(vars):
                
                pos = pos_start + 4*sin + var_indx[i]*self.np*sfl
                
                file.seek( pos )
                
                buff = read_flt_bin( file.read(self.np*sfl), sfl )

                self.df[var] = buff
                
                  
            
        # calculate the block data size in byte
        
        self.size = self.np*self.n_var*sfl + 4*sin
        
        # move file pointer to the end of current block

        file.seek( pos_start + self.size)
        