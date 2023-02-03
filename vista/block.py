#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   block.py
@Time    :   2023/02/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Class block
'''

import sys

import numpy             as     np

sys.path.append('..')

from utils.read_binary   import read_int_bin

from utils.read_binary   import read_flt_bin

from utils.read_binary   import read_log_bin


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
# - initialize an object of block and read in block data
#
# ----------------------------------------------------------------------

class BlockData:
    
    def __init__( self, file, n_var ):
        
        self.num = 0
        self.size = 0
        self.n_var = n_var
        self.dim = []
        
        # empty list for future use
        self.mean = []
        self.cor2 = []
        self.cc_vf = []
        
        # empty BlockGrid
        
        self.g = BlockGrid()

        # size of int 
        
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
        
        # read primitive variables + 'p T mu'
        
        tmp = read_flt_bin( file.read(self.np*8*sfl), sfl)
        
        self.mean = np.reshape( tmp, (8,npz,npy,npx) ).T
        
        print(type(self.mean))
        print('self.mean[:,0,0,0]')
        print(self.mean[:,0,0,0])
        print('self.mean[0,:,0,0]')
        print(self.mean[0,:,0,0])
        print('self.mean[0,0,:,0]')
        print(self.mean[0,0,:,0])
        print('self.mean[0,0,0,:]')
        print(self.mean[0,0,0,:])
        
        
        # read double correlations 
        
        tmp = read_flt_bin( file.read(self.np*36*sfl), sfl )
        
        self.cor2 = np.reshape( tmp, (36,npz,npy,npx) ).T
        
        # calculate the block data size in byte
        
        self.size = self.np*self.n_var*sfl + 4*sin

        

        