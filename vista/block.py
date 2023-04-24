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
# - initialize an object of block data.
# - for selected block, block data is read and RS is calculated.
#
# ----------------------------------------------------------------------

class BlockData:
    
    def __init__( self, file, n_var, fill ):
        
        self.verbose = False
        
        # index of block, can be read from blockdata itself
        self.num = 0
        
        # size of this BlockData (in bytes)
        self.size = 0
        self.n_var = n_var
        self.dim = []
        
        # empty list for future use
        self.mean = []
        self.cor2 = []
        self.cc_vf = []

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
        if self.num in fill: self.to_fill = True

        # primitive variables 'rho, rhou, rhov, rhow, rhoE'
        # read primitive variables + 'p T mu', when fill is true.
        if self.to_fill:
            
            tmp = read_flt_bin( file.read(self.np*8*sfl), sfl)
            
            self.mean = np.reshape( tmp, (8,npz,npy,npx) ).T
            
            if self.verbose:
                print(type(self.mean))
                print('self.mean[:,0,0,0]')
                print(self.mean[:,0,0,0])
                print('self.mean[0,:,0,0]')
                print(self.mean[0,:,0,0])
                print('self.mean[0,0,:,0]')
                print(self.mean[0,0,:,0])
                print('self.mean[0,0,0,:]')
                print(self.mean[0,0,0,:])
            
            # read double correlations: #36 is the correlations' number 
            
            tmp = read_flt_bin( file.read(self.np*36*sfl), sfl )
            
            self.cor2 = np.reshape( tmp, (36,npz,npy,npx) ).T
            
            # calculate Reynolds Stress 
            # uu represent <u`u`> actually
            
            u = self.mean[:,:,:,0]
            v = self.mean[:,:,:,1]
            w = self.mean[:,:,:,2]            
            
            self.uu = self.cor2[:,:,:,0] - np.multiply(u,u)
            self.uv = self.cor2[:,:,:,1] - np.multiply(u,v)
            self.vv = self.cor2[:,:,:,8] - np.multiply(v,v)
            self.ww = self.cor2[:,:,:,15]- np.multiply(w,w)         
        
            print("Block %d data is read."%self.num)
        # skip data chunk if fill is False
        else:
            
            # skip mean data
            file.seek( self.np*8*sfl, 1 )
            
            # skip correlation data
            file.seek( self.np*36*sfl, 1 )
        
        # calculate the block data size in byte
        
        self.size = self.np*self.n_var*sfl + 4*sin
        