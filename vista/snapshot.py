# -*- coding: utf-8 -*-
'''
@File    :   snapshot.py
@Time    :   2023/04/23 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Class and methods of reading snapshots.
             Adapted from Luis Laguarda's RAPTOR
             Further work: select only fluid cells to do dmd
'''

import os

import sys

import numpy             as     np

import pandas            as     pd

import matplotlib.pyplot as     plt

from   .read_binary      import read_int_bin

from   .read_binary      import read_flt_bin

from   .read_binary      import read_log_bin

from   .read_binary      import read_bin_3Dflt

from   .timer            import timer

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
        
        # List of position pointers to var_data chunk start
        
        self.pos_var_start = []
  
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
        
        # Name list string of variables
        
        self.vars_name = ''

        # Number of blocks
        
        self.n_bl = 0
        
        # List of blocks numbers
        
        self.bl_nums = []
        
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
#   - reading one snapshot's header, block data
#   - dropping ghost cells and assembling clean block data into 
#     a snapshot
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


        # Inform user

        if self.verbose: print( '%d blocks were read.'%self.n_bl )
        
        if self.verbose: print( 'Completed.' ); sys.stdout.flush()
        
        
        print("snapshot %d is read."%self.itstep)
        


# ----------------------------------------------------------------------
# >>> Read snapshot structure                                  ( 0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/06  - created
#
# Desc
#
#   - reading one snapshot's header, body 
#   - write snapshot structure into file for parallel reading
#   - write snapshot info for parallel dmd
#
# ----------------------------------------------------------------------

    def get_snapshot_struct( self ):

       
        # Read snapshot file

        self.read_snapshot()  
        

        # Get  the list of '| bl_num | pos_var_start | datachunk_size |' 
        
        # the number of floats, numbers of cells in x,y,z with ghost cells
        
        size = []       
        N1   = []     
        N2   = []
        N3   = []
        
        for snap_bl in self.snap_data:
            
            # size = N1 * N2 * N3 * n_var * self.kind
            
            size.append( snap_bl[1] * snap_bl[2] * snap_bl[3]
                        * self.n_vars * self.kind )
            
            N1.append( snap_bl[1] )
            N2.append( snap_bl[2] )
            N3.append( snap_bl[3] )
        
        
        # Output | bl_num | pos_var_start | datachunk_size | N1 | N2 | N3 |
        
        struct_data = np.stack(( self.bl_nums, 
                                 self.pos_var_start, 
                                 size,
                                 N1,
                                 N2,
                                 N3                   ))
        
        df_header = ['bl_num','pos','size','N1','N2','N3']
        
        df = pd.DataFrame( struct_data.T, columns=df_header )
        
        df.to_csv( 'snap_struct.csv', sep=' ',index=False )     
        
        print( f'snap_struct of snapshot {self.itstep} is output.' )
        
        
        # Write snapshots info into a file
        
        with open('snap_info.dat','w') as fi:
            
            fi.write(f'kind      {self.kind}\n')
            
            fi.write(f'n_bl      {self.n_bl}\n')
            
            fi.write(f'snap_type {self.type}\n')
            
            if self.type == 'slice':
                
                fi.write(f'slic_type {self.slic_type}\n')
            
            fi.write(f'vars_name {self.vars_name}\n')
            

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
  
  
#        if not self.snap_lean   : self.snap_lean     = True
#  
#        if not self.compressible: self.compressible  = True
#  
#        if not self.snap_with_gx: self.snap_with_gx  = True
#  
#        if not self.snap_with_tp: self.snap_with_tp  = True

  
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
  
  
        
        # Count variables
  
        if self.snap_lean: 
            
            self.n_vars += 3
            
            self.vars_name += 'u v w '
  
        else: 
            
            self.n_vars += 5
        
            self.vars_name += 'u v w xx xx '
            
  
        if self.snap_with_tp: 
            
            self.n_vars += 2
            
            self.vars_name += 'T p '
            
  
        if self.snap_with_vp: 
            
            self.n_vars += 1
        
            self.vars_name += 'vapor '
            
            
        if self.snap_with_cp: 
            
            self.n_vars += 2
            
            self.vars_name += 'cappa cp '
            
  
        if self.snap_with_mu: 
            
            self.n_vars += 1
        
            self.vars_name += 'mu '
        
        
        # Special case - skin-friction
  
        if self.snap_with_cf and self.type == 'slice':
  
            if self.slic_type == 'W': 
                
                self.n_vars += 1
                
                self.vars_name += 'cf '
            
            
        if self.snap_with_wd: 
            
            self.n_vars += 1
            
            self.vars_name += 'wd '
            

  
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
            
            print( ' - Cf     := %d' %(self.snap_with_cf) )
            
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
#     - data are stored in self.snap_data list
#     - [ bl_num, N1, N2, N3, GX, sol[n_var,:] ]
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
        
        self.bl_nums.append( bl_num )
        
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
            # if read in 3D blocks
    
            if self.snap_with_gx and self.type=='block':
                
                              
                # First direction
    
                G1  = read_bin_3Dflt( self.pos, file, N1, 1, 1, self.kind )

                self.pos += N1*self.kind

                GX.append( G1 )
                
    
                # Second direction
    
                G2  = read_bin_3Dflt( self.pos, file, 1, N2, 1, self.kind )
    
                self.pos += N2*self.kind
    
                GX.append( G2 )
                
    
                # Third direction
    
                G3  = read_bin_3Dflt( self.pos, file, 1, 1, N3, self.kind )

                self.pos     += N3*self.kind

                GX.append( G3 )
                

                    # Sanity check
    
#                    if np.all( ( G3 == 0.0 ) ): io_alert = True 
 

                # Sanity check
 
#                if np.all( ( G1 == 0.0 ) ) \
#                    or np.all( ( G2 == 0.0 ) ): io_alert = True

            # if read in slice snapshot
            
            elif self.snap_with_gx and self.type == 'slice':
                
                # slice normal to X
                
                if self.slic_type == 'X':
                    
                    
                    G2  = read_bin_3Dflt( self.pos, file, 1, N2, 1, self.kind )
                    
                    self.pos += N2*self.kind
                    
                    GX.append( G2 )
                    
                    
                    G3  = read_bin_3Dflt( self.pos, file, 1, 1, N3, self.kind )
                    
                    self.pos += N3*self.kind
                    
                    GX.append( G3 )
                    
                    
                # slice normal to Y
                
                elif self.slic_type == 'Y' or self.slic_type == 'W':
                    
                    G1  = read_bin_3Dflt( self.pos, file, N1, 1, 1, self.kind )
                    
                    self.pos += N1*self.kind
                    
                    GX.append( G1 )
                    

                    G3  = read_bin_3Dflt( self.pos, file, 1, 1, N3, self.kind )
                    
                    self.pos += N3*self.kind
                    
                    GX.append( G3 )
                
                
                # slice normal to Z
                
                elif self.slic_type == 'Z': 

                    G1  = read_bin_3Dflt( self.pos, file, N1, 1, 1, self.kind )
                    
                    self.pos += N1*self.kind
                    
                    GX.append( G1 )
                    
                    
                    G2  = read_bin_3Dflt( self.pos, file, 1, N2, 1, self.kind )
                    
                    self.pos += N2*self.kind
                    
                    GX.append( G2 )
 
 
 
            # Initialize solution buffer
 
            sol = np.zeros( shape=(self.n_vars, N1*N2*N3), dtype=np.float32 )
 
            
            # Record the starting position of var chunk
            
            self.pos_var_start.append( self.pos )
            
            
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
# >>> Drop ghost point data                                     ( 4 )
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
#
# ----------------------------------------------------------------------
   
    def drop_ghost( self , buff ):

        # Check if data is available
        
        if len(self.snap_data) == 0:
            
            raise  ValueError('please read in snapshot data first!')
        
        # When data is loaded
            
        else:    
            
            # Clean data(with out ghost cells)
            
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
                
                if self.type == 'block':
                    
                    Nx     = N1 - buff*2
                    Ny     = N2 - buff*2
                    Nz     = N3 - buff*2
                    
                
                    # Grid buffer
                    
                    GX = []
                    
                    # Loop over gx(:), gy(:), gz(:)
                    
                    for Gi in snap_data_bl[4]: GX.append( Gi[buff:-buff] )
                        
                    
                    # Solution buffer to get solution field without ghost cells
                    
                    # Notice the binary data storage order [n_var, Z, Y, X]
                    # Reshape to chunk -> slice -> reshape to vector

                    sol_buff = np.array(snap_data_bl[5]).\
                                  reshape( self.n_vars, N3, N2, N1 )
                    
                    sol_buff = sol_buff[ :, buff:-buff, buff:-buff, buff:-buff ]
                
                
                elif self.type == 'slice':
                    
                    if self.slic_type == 'X':
                        
                        Nx     = 1
                        Ny     = N2 - buff*2
                        Nz     = N3 - buff*2
                        
                        # Grid buffer
                        
                        GX = []
                        
                        # Loop over gy(:), gz(:)
                        
                        for Gi in snap_data_bl[4]: GX.append( Gi[buff:-buff] )

                        sol_buff = np.array(snap_data_bl[5]).\
                                    reshape( self.n_vars, N3, N2 )
                        
                        sol_buff = sol_buff[ :, buff:-buff, buff:-buff ]

                    
                    if self.slic_type == 'Y' or self.slic_type == 'W':
                        
                        Nx     = N1 - buff*2
                        Ny     = 1
                        Nz     = N3 - buff*2
                        
                        # Grid buffer
                        
                        GX = []
                        
                        # Loop over gx(:), gz(:)
                        
                        for Gi in snap_data_bl[4]: GX.append( Gi[buff:-buff] )

                        sol_buff = np.array(snap_data_bl[5]).\
                                    reshape( self.n_vars, N3, N1 )
                        
                        sol_buff = sol_buff[ :, buff:-buff, buff:-buff ]
                    
                    
                    if self.slic_type == 'Z':
                        
                        Nx     = N1 - buff*2
                        Ny     = N2 - buff*2
                        Nz     = 1
                        
                        # Grid buffer
                        
                        GX = []
                        
                        # Loop over gx(:), gy(:)
                        
                        for Gi in snap_data_bl[4]: GX.append( Gi[buff:-buff] )

                        sol_buff = np.array(snap_data_bl[5]).\
                                    reshape( self.n_vars, N3, N2 )
                        
                        sol_buff = sol_buff[ :, buff:-buff, buff:-buff ]
                    
                
                # Reshape the matrix ( no matter 2D or 3D) to long vectors
                       
                sol_buff = sol_buff.reshape(( self.n_vars, Nx*Ny*Nz ))
                
                # Append to snap_data_clean
                
                self.snap_cleandata.append([bl_num, Nx, Ny, Nz, GX, sol_buff])  
                        
            
            # Release memory occupied by snap_data
            
            del self.snap_data
            
            
            # Count total number of cells in the snapshots
            
            self.n_cells = 0
            
            for snap_bl in self.snap_cleandata:
                
                self.n_cells += snap_bl[1] * snap_bl[2] * snap_bl[3]

# ----------------------------------------------------------------------
# >>> check the range of snapshot                                ( 5 )
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
#
# ----------------------------------------------------------------------

    def check_range( self, clean=True ):
        
        # Check if data is available
    
            
        # Data without ghost cells is clean data
        
        if clean:
            
            if len(self.snap_cleandata) == 0:
                
                raise ValueError('please clean data first!')
            
            else:
                
                all_block_data = self.snap_cleandata
        
        else: 
            
            if len(self.snap_data) == 0:
                
                raise ValueError('please read in data first')

            else:
                
                all_block_data = self.snap_data
        
        
        # Set original range
        
        xmin =  999999.0;   xmax = -999999.0
        ymin =  999999.0;   ymax = -999999.0
        zmin =  999999.0;   zmax = -999999.0
        
        
        # Check range based on different snapshot type
        
        if self.type == 'block':
        
            for bl_data in all_block_data:
                
                xmin = min( min(bl_data[4][0]), xmin )
                xmax = max( max(bl_data[4][0]), xmax )
                ymin = min( min(bl_data[4][1]), ymin )
                ymax = max( max(bl_data[4][1]), ymax )
                zmin = min( min(bl_data[4][2]), zmin )
                zmax = max( max(bl_data[4][2]), zmax )
                
            print('snapshot range in three dimensions:')
            
            print( 'x range [ %f, %f ]' % (xmin,xmax) )
            print( 'y range [ %f, %f ]' % (ymin,ymax) )
            print( 'z range [ %f, %f ]' % (zmin,zmax) )
        
        
        elif self.type == 'slice':
            
            if self.slic_type == 'X':
                
                for bl_data in all_block_data:
                    
                    ymin = min( min(bl_data[4][0]), ymin )
                    ymax = max( max(bl_data[4][0]), ymax )
                    zmin = min( min(bl_data[4][1]), zmin )
                    zmax = max( max(bl_data[4][1]), zmax )
                    
                print('snapshot range in two dimensions:')
                
                print( 'y range [ %f, %f ]' % (ymin,ymax) )
                print( 'z range [ %f, %f ]' % (zmin,zmax) )
            
            
            elif self.slic_type == 'Y' or self.slic_type == 'W':
                
                for bl_data in all_block_data:
                    
                    xmin = min( min(bl_data[4][0]), xmin )
                    xmax = max( max(bl_data[4][0]), xmax )
                    zmin = min( min(bl_data[4][1]), zmin )
                    zmax = max( max(bl_data[4][1]), zmax )
                    
                print('snapshot range in two dimensions:')
                
                print( 'x range [ %f, %f ]' % (xmin,xmax) )
                print( 'z range [ %f, %f ]' % (zmin,zmax) )
            
            elif self.slic_type == 'Z':

                for bl_data in all_block_data:
                    
                    xmin = min( min(bl_data[4][0]), xmin )
                    xmax = max( max(bl_data[4][0]), xmax )
                    ymin = min( min(bl_data[4][1]), ymin )
                    ymax = max( max(bl_data[4][1]), ymax )
                    
                print('snapshot range in two dimensions:')
                
                print( 'x range [ %f, %f ]' % (xmin,xmax) )
                print( 'y range [ %f, %f ]' % (ymin,ymax) )



# ----------------------------------------------------------------------
# >>>  Assemble block snap_cleandata                              ( 6 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/28  - created
#
# Desc
#   - assembling cleandata(without ghost cells) blocks into a
#     snapshot's full data, return a pandas data frame.
#
# ----------------------------------------------------------------------

    def assemble_block( self ):
        
        # Different ways of assemble blocks with different shapes
        
        if self.type == 'block':
            
            bl_number = []
            x = []
            y = []
            z = []
            
            
            # Compose long vectors of coordinates x,y,z
            
            for snap_bl in self.snap_cleandata:
                
                
                bl_number.append( snap_bl[0] )
                
                x_bl = snap_bl[4][0]
                y_bl = snap_bl[4][1]
                z_bl = snap_bl[4][2]
                
                # Notice the order of output X,Y,Z !
                # https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html
                
                X,Y,Z = np.meshgrid( x_bl, y_bl, z_bl, indexing='ij' )
                
                # x,y,z are lists of nparrays
                
                x.append( np.ravel( X.T ) )
                y.append( np.ravel( Y.T ) )
                z.append( np.ravel( Z.T ) )
                
            x = np.array( x ).ravel()
            y = np.array( y ).ravel()
            z = np.array( z ).ravel()
            
            GX = np.stack( [ x.ravel(), y.ravel(), z.ravel() ] )
            
            GX_header = 'x y z '
            
            
            # Compose long vectors of solutions
            
            sol_bl = np.zeros( (self.n_vars,self.n_cells), dtype=np.float32 )
        
            pos_s = 0
            
            for snap_bl in self.snap_cleandata:
                
                pos_e = pos_s + snap_bl[1]*snap_bl[2]*snap_bl[3]
                
                sol_bl[:,pos_s:pos_e] = np.array( snap_bl[5] )
                
                pos_s = pos_e
                
                
                
        elif self.type == 'slice':
            
            
            if self.slic_type == 'X':
                
                bl_number = []
                y = []
                z = []
                
                
                # Compose long vectors of coordinates y,z
                
                for snap_bl in self.snap_cleandata:
                    
                    
                    bl_number.append(snap_bl[0])
                    
                    # GX => snap_bl[4] 
                    # For slice, only two coordinates vectors
                    
                    y_bl = snap_bl[4][0]
                    z_bl = snap_bl[4][1]
                    
                    Y, Z = np.meshgrid( y_bl, z_bl )
                    
                    y.append( np.ravel( Y ) )
                    z.append( np.ravel( Z ) )

                y = np.array( y ).ravel()
                z = np.array( z ).ravel()
                
                GX = np.stack( [ y.ravel(), z.ravel() ] )
                
                GX_header = 'y z '
                

                # Compose long vectors of solutions
                
                sol_bl = np.zeros( (self.n_vars,self.n_cells), dtype=np.float32 )
                
                pos_s = 0 
                
                for snap_bl in self.snap_cleandata:
                                        
                    pos_e = pos_s + snap_bl[1]*snap_bl[2]*snap_bl[3]
                        
                    sol_bl[:,pos_s:pos_e] = np.array( snap_bl[5] )
                        
                    pos_s = pos_e
                    
                    

            elif self.slic_type == 'W' or self.slic_type == 'Y':
                
                bl_number = []
                x = []
                z = []
                
                
                # Compose long vectors of coordinates x,z
                
                for snap_bl in self.snap_cleandata:
                    
                    
                    bl_number.append(snap_bl[0])
                    
                    x_bl = snap_bl[4][0]
                    z_bl = snap_bl[4][1]
                    
                    X, Z = np.meshgrid( x_bl, z_bl )
                    
                    x.append( np.ravel( X ) )
                    z.append( np.ravel( Z ) )

                x = np.array( x ).ravel()
                z = np.array( z ).ravel()
                
                GX = np.stack( [ x.ravel(), z.ravel() ] )
                
                GX_header = 'x z '
                

                # Compose long vectors of solutions
                
                sol_bl = np.zeros( (self.n_vars,self.n_cells), dtype=np.float32 )
                
                pos_s = 0 
                
                for snap_bl in self.snap_cleandata:
                                        
                    pos_e = pos_s + snap_bl[1]*snap_bl[2]*snap_bl[3]
                        
                    sol_bl[:,pos_s:pos_e] = np.array( snap_bl[5] )
                        
                    pos_s = pos_e


            
            elif self.slic_type == 'Z':
                
                bl_number = []
                x = []
                y = []
                
                
                # Compose long vectors of coordinates x,y
                
                for snap_bl in self.snap_cleandata:
                    
                    
                    bl_number.append(snap_bl[0])
                    
                    x_bl = snap_bl[4][0]
                    y_bl = snap_bl[4][1]
                    
                    X, Y = np.meshgrid( x_bl, y_bl )
                    
                    x.append( np.ravel( X ) )
                    y.append( np.ravel( Y ) )

                x = np.array( x ).ravel()
                y = np.array( y ).ravel()
                
                GX = np.stack( [ x.ravel(), y.ravel() ] )
                
                GX_header = 'x y '
                

                # Compose long vectors of solutions
                
                sol_bl = np.zeros( (self.n_vars,self.n_cells), dtype=np.float32 )
                
                pos_s = 0 
                
                for snap_bl in self.snap_cleandata:
                                        
                    pos_e = pos_s + snap_bl[1]*snap_bl[2]*snap_bl[3]
                        
                    sol_bl[:,pos_s:pos_e] = np.array( snap_bl[5] )
                        
                    pos_s = pos_e
            
        
        # Return pandas dataframe
        
        
        df = pd.DataFrame( GX.T, columns=GX_header.strip().split() )
        
        df_sol = pd.DataFrame( sol_bl.T, columns=self.vars_name.strip().split() )
        
        df = pd.concat([df, df_sol], axis=1)
        
        
        # Sort data frame based on coordinate z,y,x (x changes fastest)
        
        if self.type == 'block': 
            
            df.sort_values(by=['z','y','x'],inplace=True)
        
        
        elif self.type == 'slice':
            
            if self.slic_type == 'X':
                
                df.sort_values(by=['z','y'],inplace=True)
                
            
            elif self.slic_type == 'Y' or self.slic_type == 'W':
                
                df.sort_values(by=['z','x'],inplace=True)
                
            
            elif self.slic_type == 'Z':
                
                df.sort_values(by=['y','x'],inplace=True)
        
        
        self.df = df

   
                  
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
        
    test_file = test_dir1 + '/snapshot_W_002.bin'
    
    os.chdir( test_dir1 )
    
    snapshot1 = Snapshot( test_file )
    
#    snapshot1.read_snapshot()
    
    with timer('get snapshot struct'):
        
        snapshot1.get_snapshot_struct()
        
    
#    snapshot1.verbose = False

'''
    with timer('read one snapshot '):
    
        snapshot1.read_snapshot()

        # Drop ghost points 
        
        snapshot1.drop_ghost( buff=3 )
        
        
        # Assemble blocks data to one whole snapshot
        # (sorting is included in self.assemble_block() )
        
        snapshot1.assemble_block()
        
        
        print(snapshot1.df)
        
    with timer('show one slice '):
        
        x = np.array(snapshot1.df['x'])
        
        y = np.array(snapshot1.df['y'])
        
        z = np.array(snapshot1.df['z'])
        
        p = np.array(snapshot1.df['p'])
        
#        print(np.unique(x))
        
        N_x = len(np.unique(x))
        N_y = len(np.unique(y))
        N_z = len(np.unique(z))
        
        x = x.reshape(N_z,N_y,N_x)
        y = y.reshape(N_z,N_y,N_x)
        z = z.reshape(N_z,N_y,N_x)
        p = p.reshape(N_z,N_y,N_x)
        
        x = x[50,:,:]
        y = y[50,:,:]
        p = p[50,:,:]
        
        
        fig, ax = plt.subplots()
        contour = ax.pcolor(x,y,p)
        ax.set_title('pressure')
        plt.show()
'''
'''   
    with timer('check one data block:'):
        
        # [bl_num, Nx, Ny, Nz, GX, sol_buff]
        
        snap_bl = snapshot1.snap_cleandata[0]
        
        Nx = snap_bl[1]
        Ny = snap_bl[2]
        Nz = snap_bl[3]
        
        print(Nx,Ny,Nz)
        
        x = snap_bl[4][0]
        z = snap_bl[4][1]
        X,Z = np.meshgrid( x, z )
        
        print(len(x),len(z))
        
        u = np.array(snap_bl[5][0]).reshape(Nz,Nx)
 
        u_slice = u.T     
        
        print(X.shape)
        print(Z.shape)
        print(u_slice.shape)   
        
        fig, ax = plt.subplots()
        contour = ax.contourf(X.T,Z.T, u_slice)
        ax.set_title('Contour Plot')
        plt.show()
'''         
        
'''        

        
#        for l in u_slice:  print(l)

        X,Y = np.meshgrid( x, y )

        mesh = pv.StructuredGrid(X,Y,u_slice).elevation()
        
        contour = mesh.contour()
        
        plotter = pv.Plotter()
        
        plotter.add_mesh(contour, line_width=1, color='white')
        
        plotter.show()

        
        print(type(X),type(Y))
        print(X.shape,Y.shape)
        
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
