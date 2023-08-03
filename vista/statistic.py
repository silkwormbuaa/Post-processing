# -*- coding: utf-8 -*-
'''
@File    :   statistic.py
@Time    :   2022/10/12 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os

from   .block            import BlockData

from   .io_binary        import read_int_bin

from   .io_binary        import read_flt_bin

from   .io_binary        import read_log_bin

#%% 
class StatisticData:

# ----------------------------------------------------------------------
# >>> INITIALIZE CLASS INSTANCE                                   ( 0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/12  - created
#
# Input
#     
# - statistic binary file directory
# ----------------------------------------------------------------------

    def __init__( self, stat_dir ):
        
        # directory to the statistic file
        
        self.dir = stat_dir

        # file size
        
        self.fsize = os.stat(self.dir).st_size
        
        # file pointer position
        
        self.pos = 0

        # number of variables
        
        self.n_var = 0
        
        # list of blocks
        
        self.bl = []
               
        # Verbose ? 
        
        self.verbose = False
#%%
# ----------------------------------------------------------------------
# >>> READ STATISTIC BINARY FILE HEADER                           ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/12  - created
#
# Input
# 
# - opened statistic file
# - 
# ----------------------------------------------------------------------

    def read_stat_header( self, file ):

        # Internal variables

        sfl = 8
        sin = 4
        slg = 4

        # Read file format

        self.format = read_int_bin(file.read(sin),sin)

        self.pos += sin


        if self.format == 1:

            # Read number of samples

            self.n_samples = read_int_bin(file.read(sin),sin)

            self.pos += sin

            # Read sample step & start step

            self.sample_step = read_int_bin(file.read(sin),sin)
            self.start_step  = read_int_bin(file.read(sin),sin)

            self.pos += int(2*sin)

            # Read sample time & start time

            self.sample_time = read_flt_bin(file.read(sfl),sfl)
            self.start_time  = read_flt_bin(file.read(sfl),sfl)
            
            self.pos += int(2*sfl)

            # Read settings; slg - size of logical

            buf_log = read_log_bin(file.read(12*slg),slg)

            self.pos += int(12*slg)

            # Header size

            self.header_size = self.pos
            
            # Store settings
            
            self.meanvalues       = buf_log[ 0]
            if self.meanvalues: self.n_var += 8

            self.doublecorr       = buf_log[ 1]
            if self.doublecorr: self.n_var += 36

            self.triplecorr       = buf_log[ 2]

            self.quadruplecorr    = buf_log[ 3]

            self.autocorr         = buf_log[ 4]

            self.mean_invar       = buf_log[ 5]

            self.schlieren        = buf_log[ 6]

            self.cavitation_stats = buf_log[ 7]

            self.vapor_gas_stats  = buf_log[ 8]

            self.rste             = buf_log[ 9]

            self.thermo           = buf_log[10]

            self.visc_diff        = buf_log[11]
            
            # if want to see header, set verbose True.            
            if (self.verbose is True) : 
                print('format '          + str(self.format)           + '\n')
                print('n_samples '       + str(self.n_samples)        + '\n')
                print('start_step '      + str(self.start_step)       + '\n')
                print('sample_step '     + str(self.sample_step)      + '\n')
                print('sample_time '     + str(self.sample_time)      + '\n')
                print('start_time '      + str(self.start_time)       + '\n')
                print('meanvalues '      + str(self.meanvalues)       + '\n')
                print('doublecorr '      + str(self.doublecorr)       + '\n')
                print('triplecorr '      + str(self.triplecorr)       + '\n')
                print('quadruplecorr '   + str(self.quadruplecorr)    + '\n')
                print('autocorr '        + str(self.autocorr)         + '\n')
                print('mean_invar '      + str(self.mean_invar)       + '\n')
                print('schlieren '       + str(self.schlieren)        + '\n')
                print('cavitation_stats '+ str(self.cavitation_stats) + '\n')
                print('vapor_gas_stats ' + str(self.vapor_gas_stats)  + '\n')
                print('rste '            + str(self.rste)             + '\n')
                print('thermo '          + str(self.thermo)           + '\n')
                print('visc_diff '       + str(self.visc_diff)        + '\n')
                print('header_size '     + str(self.header_size)      + '\n')

        else:
            raise ValueError('Header format in statistics.bin not supported')


    #    if rank == root:

    #        print( 'STATS Number of samples: %d' %(self.n_samples))

# ----------------------------------------------------------------------
# >>> READ STATISTIC BINARY FILE BODY                             ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/12  - created
#
# Input
#
# - opened file
# - a list of block numbers whose data chunk will be filled
#   (no index issue)
# ----------------------------------------------------------------------

    def read_stat_body( self, file , fill ):
        
        end_of_file = False
               
        # New position after reading headers
        
        file.seek( self.header_size )
        
        self.pos = self.header_size
        
        while not end_of_file:
        
            # read in block one by one
            # only blocks in fill list will be filled with data chunk
        
            self.bl.append( BlockData( file, self.n_var, fill ) )
            
            self.pos = self.pos + self.bl[-1].size
                                    
            if self.pos >= self.fsize: end_of_file = True

"""
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
"""