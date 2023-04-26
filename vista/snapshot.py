# -*- coding: utf-8 -*-
'''
@File    :   snapshot.py
@Time    :   2023/04/23 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Class and methods of reading snapshots.
         :   Adapted from Luis Laguarda's RAPTOR
'''

import os

import sys

import numpy             as np

import pyvista           as pv

import matplotlib.pyplot        as plt

sys.path.append('..')

from   utils.read_binary import read_int_bin

from   utils.read_binary import read_flt_bin

from   utils.read_binary import read_log_bin

from   utils.read_binary import read_bin_3Dflt

from   utils.timer       import timer

class Snapshot:

# ----------------------------------------------------------------------
# >>> Initialize snapshot class                                  ( 0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def __init__( self, snap_dir ):
        
        # Directory
  
        self.dir = snap_dir

        # Chekc file sizes

        self.fsize = os.stat( self.dir ).st_size
        
        # Snapshot type
  
        if snap_dir[-9:-8] in [ 'W', 'X', 'Y', 'Z' ]:
  
            self.type = 'slice'; self.slic_type = snap_dir[-9:-8]
  
        else:
  
            self.type = 'block'
  
  
        # Snapshot data (per block)
  
        self.snap_data = []
  
        # Position pointers
  
        self.pos = 0
  
        # Characteristics
  
        self.kind          = 4
  
        self.itstep        = 0
  
        self.itstep_check  = -1
  
        self.itime         = 0.0
  
        self.n_species     = 0
  
        self.snap_lean     = False
  
        self.compressible  = False
  
        self.snap_with_gx  = False
  
        self.snap_with_tp  = False
  
        self.snap_with_vp  = False
  
        self.snap_with_cp  = False
  
        self.snap_with_mu  = False
  
        self.snap_with_wd  = False
  
        self.snap_with_cf  = False
  
        self.snap_with_bg  = False
  
  
        # Number of variables
  
        self.n_vars = 0  

        # Number of blocks
        
        self.n_bl = 0
        # Verbose
  
        self.verbose = True


# ----------------------------------------------------------------------
# >>> Read snapshot                                              ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_snapshot( self ):
        
        # End of file marker

        end_of_file = False


# --- Read snapshot file

        with open( self.dir, 'rb' ) as fs:

            # Start from the beginning of the file

            fs.seek( self.pos )

            # Read snapshot header

            self.read_snap_header( fs )
            
            self.pos += self.header_size

            # Read body

            while not end_of_file:

                # Read current block 
                # self.pos is updated inside self.read_snap_block

                self.read_snap_block( fs )
                
                self.n_bl += 1

                # End of file?

                if self.pos >= self.fsize: end_of_file = True

        # Close file

        fs.close()

        # Inform user

        if self.verbose: print( '%d blocks were read.'%self.n_bl )
        
        if self.verbose: print( 'Completed.' ); sys.stdout.flush()

# ----------------------------------------------------------------------
# >>> Read snapshot header                                       ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_snap_header( self, file ):
  
        # Book-keeping
        
        self.header_size = 0
  
        # Default float size
  
        sin = 4
  
        slg = 4
  
        sfl = 8
  
        # Header format
  
        hformat = read_int_bin( file.read(sin), sin )
  
        # Floating point precision format
  
        kind    = read_int_bin( file.read(sin), sin )
    
        # Just in case
  
        if kind != 4 and kind != 8: kind = 4
  
        # Time step
  
        itstep  = read_int_bin( file.read(sin), sin )
  
        # Number of species
  
        nspec   = read_int_bin( file.read(sin), sin )
  
        self.header_size += 4*sin

        # Simulation time
  
        itime   = read_flt_bin( file.read(sfl), sfl )
  
        self.header_size += sfl
  
        # Content
  
        buf_log = read_log_bin( file.read(10*slg), slg )
  
        # Assign to class variables
  
        self.hformat       = hformat
  
        self.kind          = kind
  
        self.itstep        = itstep
  
        self.itime         = itime
  
        self.n_species     = nspec
  
        self.snap_lean     = buf_log[0]
  
        self.compressible  = buf_log[1]
  
        self.snap_with_gx  = buf_log[2]
  
        self.snap_with_tp  = buf_log[3]
  
        self.snap_with_vp  = buf_log[4]
  
        self.snap_with_cp  = buf_log[5]
  
        self.snap_with_mu  = buf_log[6]
  
        # just in case ?
        
#        if self.itstep != self.itstep_check:
  
#            self.itstep = self.itstep_check
  
  
        if not self.snap_lean   : self.snap_lean     = True
  
        if not self.compressible: self.compressible  = True
  
        if not self.snap_with_gx: self.snap_with_gx  = True
  
        if not self.snap_with_tp: self.snap_with_tp  = True

  
        # Block-snapshots
  
        if self.type == 'block':
  
            self.snap_with_wd  = buf_log[7]
  
            self.header_size += 8*slg
  
        # Slice-snapshots
  
        else:
  
            self.snap_with_cf  = buf_log[7]
  
            self.snap_with_wd  = buf_log[8]
  
            self.snap_with_bg  = buf_log[9]
  
            self.header_size += 10*slg
  

            if self.slic_type == 'W' and not self.snap_with_cf:
  
                self.snap_with_cf = True
  
            if self.slic_type == 'Z' and not self.snap_with_wd:
  
                self.snap_with_wd = True
  
        # Count variables
  
        if self.snap_lean   : self.n_vars += 3
  
        else                : self.n_vars += 5
  
        if self.snap_with_tp: self.n_vars += 2
  
        if self.snap_with_wd: self.n_vars += 1
  
        if self.snap_with_vp: self.n_vars += 1
  
        if self.snap_with_cp: self.n_vars += 2
  
        if self.snap_with_mu: self.n_vars += 1
  

        # Special case - skin-friction
  
        if self.snap_with_cf and self.type == 'slice':
  
            if self.slic_type == 'W': self.n_vars += 1
  
  
        # Inform user
  
        if self.verbose:
  
            print( '' )
   
            print( 'Current snapshot: ...' + self.dir[-50:] )
   
            print( ' - File size is %d Mb (%d B)'%(self.fsize/(1000000),self.fsize) )
   
            print( ' - Kind   := %d' %(self.kind  ) )
   
            print( ' - Fields := %d' %(self.n_vars) )
   
            print( ' - Time   := %.4e'%(self.itime ) )
   
            print( ' - Step   := %d' %(self.itstep) )
            
            print( ' - lean   := %d' %(self.snap_lean) )
            
            print( ' - compressible := %d' %(self.compressible) )
            
            print( ' - gx     := %d' %(self.snap_with_gx) )
            
            print( ' - tp     := %d' %(self.snap_with_tp) )
   
            print( ' - vp     := %d' %(self.snap_with_vp) )
            
            print( ' - cp     := %d' %(self.snap_with_cp) )
            
            print( ' - mu     := %d' %(self.snap_with_mu) )
            
            print( ' - wd     := %d' %(self.snap_with_wd) )
            
            print( ' - n_vars := %d' %(self.n_vars) )
            
            print( '' ); sys.stdout.flush()

# ----------------------------------------------------------------------
# >>> Read snapshot block                                        ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/23  - created
#
# Desc
#     - read one snapshot block data
# ----------------------------------------------------------------------

    def read_snap_block( self, file ):

        # Size of int

        sin = 4
 
        # Verbose
 
        verbose = self.verbose
 
        # New position in file
 
        file.seek( self.pos )
 
        # Read global block number and block dimensions
 
        bl_num = read_int_bin( file.read(   sin ), sin )
 
        bl_dim = read_int_bin( file.read( 3*sin ), sin )
 
        # This is a sanity check - related to an old IO bug
 
        if bl_num == 0:
 
            # Block dimensions
 
            N1 = self.snap_data[-1][1]
 
            N2 = self.snap_data[-1][2]
 
            N3 = self.snap_data[-1][3]
 
            # Use previous results
 
            N3n = N3
 
            if N3 == 1: N3n = 0
 
            # Update block size
 
            bl_size = N1 + N2 + N3n + (N1*N2*N3)*self.n_vars
 
            # Update file pointer
 
            self.pos += 4*sin + bl_size * self.kind
 
            # Inform user
 
            if verbose:
 
                str_ = ' Faulty block %5d - %2d x %2d x %2d'%(bl_num,N1,N2,N3)
 
                print( str_ )
 
                sys.stdout.flush()
 
        else:
                
            # Corrupted IO flag
 
            io_alert = False
 
            # Block dimensions
 
            N1 = bl_dim[0]
 
            N2 = bl_dim[1]
 
            N3 = bl_dim[2]
 
            # Update file pointer
 
            self.pos += 4*sin
 
            # Grid buffer
 
            GX = []
 
            # Read grid vectors
    
            if self.snap_with_gx:
    
                # First direction
    
                G1  = read_bin_3Dflt( self.pos, file, N1, 1, 1, self.kind )
    
                self.pos += N1*self.kind
    
                GX.append( G1 )
    
                # Second direction
    
                G2  = read_bin_3Dflt( self.pos, file, 1, N2, 1, self.kind )
    
                self.pos += N2*self.kind
    
                GX.append( G2 )
    
                # Third direction
    
                if N3 > 1:
    
                    G3  = read_bin_3Dflt( self.pos, file, 1, 1, N3, self.kind )
    
                    self.pos     += N3*self.kind
    
                    GX.append( G3 )
    
    
    
                    # Sanity check
    
                    if np.all( ( G3 == 0.0 ) ): io_alert = True 
 

                # Sanity check
 
                if np.all( ( G1 == 0.0 ) ) \
                    or np.all( ( G2 == 0.0 ) ): io_alert = True
 
            # Initialize solution buffer
 
            sol = np.zeros( shape=(self.n_vars, N1*N2*N3), dtype=np.float32 )
 
 
            # Read solution fields
 
            for v in range( self.n_vars ):
 
                sol[v,:] = read_bin_3Dflt( self.pos, file, N1,N2,N3, self.kind )
 
                self.pos += ( N1*N2*N3 ) * self.kind
 
                # Sanity check on the pressure
 
                if v == 4 and np.all( ( sol[v,:] == 0.0 ) ): io_alert = True


            # Evaluate
 
            if io_alert:

            # Inform user
 
                if verbose:
 
                    str_ = ' Faulty block %5d - %2d x %2d x %2d' \
                        %( bl_num, N1, N2, N3 ); print( str_ )
 
                    sys.stdout.flush()
 
            else:
 
                # Inform user that all checks are passed
 
                if verbose:
 
                    str_ = ' Block %5d has size %2d x %2d x %2d' \
                        %( bl_num, N1, N2, N3 ); print( str_ )
 
                    sys.stdout.flush()
 
                # Append to snap_data
 
                self.snap_data.append( [ bl_num, N1, N2, N3, GX, sol ] )


# ----------------------------------------------------------------------
# >>> check the range of snapshot                                ( 4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/24  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def check_range( self ):
        
        # check if data is available
        
        if len(self.snap_data) == 0:
            
            print('please read in snapshot data first!')
        
        # when data is loaded
            
        else:
            
            # check data with ghost cells
            
            xmin =  999999.0;   xmax = -999999.0
            ymin =  999999.0;   ymax = -999999.0
            zmin =  999999.0;   zmax = -999999.0
            
            for bl_data in self.snap_data:
                
                xmin = min( min(bl_data[4][0]), xmin )
                xmax = max( max(bl_data[4][0]), xmax )
                ymin = min( min(bl_data[4][1]), ymin )
                ymax = max( max(bl_data[4][1]), ymax )
                zmin = min( min(bl_data[4][2]), zmin )
                zmax = max( max(bl_data[4][2]), zmax )
                
            print('before dropping ghost cells')
            print( 'x range [ %f, %f ]' % (xmin,xmax) )
            print( 'y range [ %f, %f ]' % (ymin,ymax) )
            print( 'z range [ %f, %f ]' % (zmin,zmax) )

        
        # check if data is available
        
        if len(self.snap_cleandata) == 0:
            
            print('please clean snapshot data first!')
        
        # when data is loaded
            
        else:

            # check range after drop ghost cells
            
            xmin =  999999.0;   xmax = -999999.0
            ymin =  999999.0;   ymax = -999999.0
            zmin =  999999.0;   zmax = -999999.0
            
            for bl_data in self.snap_cleandata:
                
                xmin = min( min(bl_data[4][0]), xmin )
                xmax = max( max(bl_data[4][0]), xmax )
                ymin = min( min(bl_data[4][1]), ymin )
                ymax = max( max(bl_data[4][1]), ymax )
                zmin = min( min(bl_data[4][2]), zmin )
                zmax = max( max(bl_data[4][2]), zmax )
                
            print('after dropping ghost cells:')
            print( 'x range [ %f, %f ]' % (xmin,xmax) )
            print( 'y range [ %f, %f ]' % (ymin,ymax) )
            print( 'z range [ %f, %f ]' % (zmin,zmax) )


# ----------------------------------------------------------------------
# >>> Drop ghost point data                                     ( 5 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/24  - created
#
# Desc
#
# ----------------------------------------------------------------------
   
    def drop_ghost( self , buff ):

        # check if data is available
        
        if len(self.snap_data) == 0:
            
            raise  ValueError('please read in snapshot data first!')
        
        # when data is loaded
            
        else:    
            
            # clean data(with out ghost cells)
            
            self.snap_cleandata = []  
            
            
            # Loop over all snap_data's elements
            # Snap_data_bl = [ bl_num, N1, N2, N3, GX, sol ]
            
            for snap_data_bl in self.snap_data:
                
                bl_num = snap_data_bl[0]
                
                # N1, N2, N3 are dimensions with ghost cells
                
                N1     = snap_data_bl[1]
                N2     = snap_data_bl[2]
                N3     = snap_data_bl[3]
                
                # Nx, Ny, Nz are dimensions without ghost cells
                
                Nx     = N1 - buff*2
                Ny     = N2 - buff*2
                Nz     = N3 - buff*2
                
                # Grid buffer
                
                GX = []
                
                # Loop over gx(:), gy(:), gz(:)
                
                for Gi in snap_data_bl[4]: GX.append( Gi[buff:-buff] )
                    
                
                # Solution buffer to get solution field without ghost cells
                
                # Notice the binary data storage order [n_var, Z, Y, X]
                # reshape to chunk -> slice -> reshape to vector

                sol_buff = np.array(snap_data_bl[5]).reshape( self.n_vars, N3, N2, N1 )
                
                sol_buff = sol_buff[ :, buff:-buff, buff:-buff, buff:-buff ]
                
                sol_buff = sol_buff.reshape(( self.n_vars, Nx*Ny*Nz ))
                
                
                # Append to snap_data_clean
                
                self.snap_cleandata.append( [bl_num, Nx, Ny, Nz, GX, sol_buff] )
                    
                    
# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/25  - created
#
# Desc
#
# ----------------------------------------------------------------------
def Testing():
    
    test_dir1 = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots/snapshot_00514923'
    
    test_dir2 = '/home/wencanwu/my_simulation/temp/220926_lowRe/snapshots/snapshot_00452401'
    
    test_dir = test_dir2 + '/snapshot.bin'
    
    snapshot1 = Snapshot( test_dir )
    
#    snapshot1.verbose = False

    with timer('read one snapshot:'):
    
        snapshot1.read_snapshot()
                
        snapshot1.drop_ghost( buff=3 )
        
        snapshot1.check_range()
        
    with timer('check one data block:'):
        
        # [bl_num, Nx, Ny, Nz, GX, sol_buff]
        
        snap_bl = snapshot1.snap_cleandata[0]
        
        Nx = snap_bl[1]
        Ny = snap_bl[2]
        Nz = snap_bl[3]
        
        x = snap_bl[4][0]
        y = snap_bl[4][1]
        
        u = np.array(snap_bl[5][0]).reshape(Nz,Ny,Nx)
        
        u_slice = u[0,:,:]
        
        u_slice = u_slice.T
        
        print(type(u_slice))
        print(type(u_slice[0][0]))
        print(u_slice.shape)
        
#        for l in u_slice:  print(l)

'''        
        X,Y = np.meshgrid( x, y )
        
        print(type(X),type(Y))
        print(X.shape,Y.shape)
        
        fig, ax = plt.subplots()
        contour = ax.contourf(X,Y, u_slice)
#        ax.clabel(contour, inline=True, fontsize=10)
        ax.set_title('Contour Plot')
        plt.show()

'''
# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -- )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":
    
    Testing()
