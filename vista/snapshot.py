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
import gc
import sys
import copy

import numpy             as     np
import pandas            as     pd
import pyvista           as     pv
import tecplot           as     tp
from   tecplot.constant  import FieldDataType
from   copy              import deepcopy

from   .io_binary        import read_int_bin
from   .io_binary        import read_flt_bin
from   .io_binary        import read_log_bin
from   .io_binary        import read_3Dflt_bin

from   .io_binary        import write_flt_bin
from   .io_binary        import write_int_bin
from   .io_binary        import write_log_bin

from   .io_vtk           import create_3d_vtkRectilinearGrid
from   .io_vtk           import add_var_vtkRectilinearGrid
from   .io_vtk           import create_multiblock_dataset
from   .io_vtk           import write_vtm_file

from   .grid             import GridData

from   .block            import SnapBlock

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

    def __init__( self, snap_dir=None ):
        
        """
        snap_dir: (optional) if None, generate a void Snapshot
        """
        
        # Directory
        
        if snap_dir is not None:
  
            self.dir = snap_dir

            # Chekc file sizes

            self.fsize = os.stat( self.dir ).st_size
            
            # Snapshot type
    
            if snap_dir[-9:-8] in [ 'W', 'X', 'Y', 'Z' ]:
    
                self.type = 'slice'; self.slic_type = snap_dir[-9:-8]
    
            else:
    
                self.type = 'block'
        
        else:
            
            self.dir = None 
            
            self.fsize = 0
            
            self.type = None
  
        # Snapshot data (per block)
  
        self.snap_data = []
        self.snap_cleandata = []
  
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
  
        self.n_var = 0  
        
        # Name list string of variables
        
        self.vars_name = []

        # Number of blocks
        
        self.n_bl = 0
        
        # List of blocks numbers (of blocks in binary files)
        
        self.bl_nums = []
        self.bl_nums_clean = []
        
        # Verbose
  
        self.verbose = False
        
        # Grid3d (including all blocks grids from inca_grid.bin file)

        self.grid3d = None
        
        # number of actually filled blocks
        
        self.filled = 0


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
#
# ----------------------------------------------------------------------

    def read_snapshot( self, block_list=None, var_read=None ):
        
        """
        block_list: (optional) if None, read all blocks.
        """
        
        # before read in snapshots, clear first.
        
        if len(self.bl_nums) != 0: self.__init__(self.dir)
        
        # End of file marker

        end_of_file = False

# ----- Read snapshot file

        with open( self.dir, 'rb' ) as fs:

            # Start from the beginning of the file
            
            fs.seek( self.pos )

# --------- Read snapshot header
            
            self.read_snap_header( fs )
            self.pos += self.header_size
            
            # different format have different length of header
            fs.seek( self.header_size ) 

            print(f"...{self.dir[-25:]} has variable {self.vars_name}")

# --------- Read body
            
            while not end_of_file:

                # Read current block 
                # self.pos is updated inside self.read_snap_block
                 
                self.snap_data.append( SnapBlock(fs, 
                                                 block_list,
                                                 var_read, 
                                                 self.n_var,
                                                 self.vars_name,
                                                 self.snap_with_gx,
                                                 self.type) )
                
                # collect the bl_nums and start positions of solution
                self.bl_nums.append( self.snap_data[-1].num )
                self.pos_var_start.append(self.snap_data[-1].pos_var_start)
                
                # append grid to snap_data
                
                if self.grid3d is not None:
                    self.snap_data[-1].g = self.grid3d.g[self.bl_nums[-1]-1]
                
                # count the number of blocks in the snapshot file
                self.n_bl += 1
                
                # count the number of filled blocks
                if self.snap_data[-1].to_fill: self.filled += 1
                
                # verbose
                
                if self.verbose: 
                    print(f"Block {self.bl_nums[-1]}, dimension:",end='')
                    print(f"{self.snap_data[-1].npx}x",end='')
                    print(f"{self.snap_data[-1].npy}x",end='')
                    print(f"{self.snap_data[-1].npz}, ",end='')
                    print(f"filled: {self.snap_data[-1].to_fill}")
                    sys.stdout.flush()
                    
                # End of file?

                self.pos = fs.tell()
                if self.pos >= self.fsize: end_of_file = True
                

# ----- Inform user

        print( f'\n{self.n_bl} blocks were read ',end='')
        print( f'with { self.filled } blocks filled.' )

        if self.type == 'block': snap_type = self.type
        if self.type == 'slice': snap_type = self.slic_type
        
        print(f"\nSnapshot {self.itstep} {snap_type} is read.\n")
        
        sys.stdout.flush()

# ----------------------------------------------------------------------
# >>> Read snapshot and output snapshot structure/info            ( 0 )
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

    def get_snapshot_struct( self, struct_file=None, info_file=None ):
        
        """
        struct_file : if not given, by default snap_struct.csv.
        info_file   : if not given, by default snap_info.dat.
        """

        # initialize output struct/info file name
       
        if struct_file is None: 
            struct_file = 'snap_struct.csv'
       
        if info_file is None: 
            info_file = 'snap_info.dat'

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
            
            size.append( snap_bl.npx * snap_bl.npy * snap_bl.npz
                        * self.n_var * self.kind )
            
            N1.append( snap_bl.npx )
            N2.append( snap_bl.npy )
            N3.append( snap_bl.npz )
        
        
        # Output | bl_num | pos_var_start | datachunk_size | N1 | N2 | N3 |
        
        struct_data = np.stack(( self.bl_nums, 
                                 self.pos_var_start, 
                                 size,
                                 N1,
                                 N2,
                                 N3                   ))
        
        df_header = ['bl_num','pos','size','N1','N2','N3']
        
        df = pd.DataFrame( struct_data.T, columns=df_header )
        
        df.to_csv( struct_file, sep=' ',index=False )     
        
        # Write snapshots info into a file
        
        with open( info_file,'w') as fi:
            
            fi.write(f'kind      {self.kind}\n')
            
            fi.write(f'n_bl      {self.n_bl}\n')
            
            fi.write(f'snap_type {self.type}\n')
            
            if self.type == 'slice':
                
                fi.write(f'slic_type {self.slic_type}\n')
            
            fi.write(f'vars_name {self.vars_name}\n')
        
        
        print(f'snap struct/info files of snapshot {self.itstep} are output.\n')
            

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
        
        """
        file : opened file object
        """
  
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
            
            self.n_var += 3
            self.vars_name += ['u', 'v', 'w']
  
        else: 
            
            self.n_var += 5
            self.vars_name += ['u', 'v', 'w', 'xx', 'xx']
            
  
        if self.snap_with_tp: 
            
            self.n_var += 2
            self.vars_name += ['T', 'p']
            
  
        if self.snap_with_vp: 
            
            self.n_var += 1
            self.vars_name += ['vapor']
            
            
        if self.snap_with_cp: 
            
            self.n_var += 2
            self.vars_name += ['cappa', 'cp']
            
  
        if self.snap_with_mu: 
            
            self.n_var += 1
            self.vars_name += ['mu']
        
        
        # Special case - skin-friction
  
        if self.snap_with_cf and self.type == 'slice':
  
            if self.slic_type == 'W': 
                
                self.n_var += 1
                self.vars_name += ['cf']
            
            
        if self.snap_with_wd: 
            
            self.n_var += 1
            self.vars_name += ['wd']
            

        # Inform user
  
        if self.verbose:
  
            print( '' )
            print( 'Current snapshot: ...' + self.dir[-50:] )
            print( ' - File size is %d Mb (%d B)'%(self.fsize/(1000000),self.fsize) )
            print( ' - Kind   := %d' %(self.kind  ) )
            print( ' - Fields := %d' %(self.n_var) )
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
            print( ' - n_var  := %d' %(self.n_var) )
   
            print( '' ); sys.stdout.flush()
          

# ----------------------------------------------------------------------
# >>> compute gradients in snapshots                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/02/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_gradients( self, block_list=None, grads=['grad_rho'] ):
        
        """
        block_list: list of blocks that are going to compute gradients;
                    if None, all blocks will compute gradients.
        grads: list of strings, choose from 
        ['grad_rho', 'laplacian', 'grad_p', 'vorticity','Q_cr','lambda2','div'],
        default is ['grad_rho']
        """
        
        if block_list is None:
            
            for snap_bl in self.snap_data:
                snap_bl.compute_gradients_block( grads )
        
        else:
            
            for num in block_list:
                snap_bl = self.snap_data[self.bl_nums.index(num)]
                snap_bl.compute_gradients_block( grads )
            
        print(f"Snapshot {self.itstep} gradients ({grads}) are computed.\n")
             

# ----------------------------------------------------------------------
# >>> Drop ghost point data                                       (Nr.)
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

    def drop_ghost( self, block_list=None, buff=3, mode='symmetry' ):

        """
        self.snap_data will be replaced by self.snap_cleandata
        mode: 'symmetry' or 'oneside'
        
        self.n_cells_clean (number of cell after dropping ghost) will be updated
        """
        
        # Check if data is available
        
        if len(self.snap_data) == 0:
            raise  ValueError('please read in snapshot data first!')
        
        # when data is loaded
        
        else:
            
            del self.snap_cleandata
            gc.collect()
            
            # Clean data(with out ghost cells)
            self.snap_cleandata = []
            
            if block_list is None:
                block_list = self.bl_nums
            
            for block in self.snap_data:
                
                if block.num not in block_list:
                    continue
                
                self.snap_cleandata.append( block.drop_ghost(buff=buff,mode=mode) )

        # update list of block numbers
        # count total number of cells in the snapshots after dropping ghost
        
        self.bl_nums_clean = [bl.num for bl in self.snap_cleandata]
        
        self.n_cells_clean = 0
        for bl in self.snap_cleandata:
            self.n_cells_clean += bl.npx * bl.npy * bl.npz
        
        print(f"Snapshot {self.itstep} dropped ghost cells.")
        

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
        
        """
        print x,y,z range of data
        """
    
            
        # Data without ghost cells is clean data
        
        if clean:
            
            if len(self.snap_cleandata) == 0:
                
                raise ValueError('please clean data first!')
            
            else:
                
                all_blocks = self.snap_cleandata
        
        else: 
            
            if len(self.snap_data) == 0:
                
                raise ValueError('please read in data first')

            else:
                
                all_blocks = self.snap_data
        
        
        # Set original range
        
        xmin =  999999.0;   xmax = -999999.0
        ymin =  999999.0;   ymax = -999999.0
        zmin =  999999.0;   zmax = -999999.0
        
        
        # Check range based on different snapshot type
        
        if self.type == 'block':
        
            for bl_data in all_blocks:
                
                xmin = min( min(bl_data.npx), xmin )
                xmax = max( max(bl_data.npx), xmax )
                ymin = min( min(bl_data.npy), ymin )
                ymax = max( max(bl_data.npy), ymax )
                zmin = min( min(bl_data.npz), zmin )
                zmax = max( max(bl_data.npz), zmax )
                
            print('snapshot range in three dimensions:')
            
            print( 'x range [ %f, %f ]' % (xmin,xmax) )
            print( 'y range [ %f, %f ]' % (ymin,ymax) )
            print( 'z range [ %f, %f ]' % (zmin,zmax) )
        
        
        elif self.type == 'slice':
            
            if self.slic_type == 'X':
                
                for bl_data in all_blocks:
                    
                    ymin = min( min(bl_data.npy), ymin )
                    ymax = max( max(bl_data.npy), ymax )
                    zmin = min( min(bl_data.npz), zmin )
                    zmax = max( max(bl_data.npz), zmax )
                    
                print('snapshot range in two dimensions:')
                
                print( 'y range [ %f, %f ]' % (ymin,ymax) )
                print( 'z range [ %f, %f ]' % (zmin,zmax) )
            
            
            elif self.slic_type == 'Y' or self.slic_type == 'W':
                
                for bl_data in all_blocks:
                    
                    xmin = min( min(bl_data.npx), xmin )
                    xmax = max( max(bl_data.npx), xmax )
                    zmin = min( min(bl_data.npz), zmin )
                    zmax = max( max(bl_data.npz), zmax )
                    
                print('snapshot range in two dimensions:')
                
                print( 'x range [ %f, %f ]' % (xmin,xmax) )
                print( 'z range [ %f, %f ]' % (zmin,zmax) )
            
            elif self.slic_type == 'Z':

                for bl_data in all_blocks:
                    
                    xmin = min( min(bl_data.npx), xmin )
                    xmax = max( max(bl_data.npx), xmax )
                    ymin = min( min(bl_data.npz), ymin )
                    ymax = max( max(bl_data.npz), ymax )
                    
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
        
        """
        self.snap_cleandata should be ready!
        assembling cleandata(without ghost cells) blocks into a
        snapshot's full data, return a pandas data frame.
        """
        
        # Different ways of assemble blocks with different shapes

        vars = self.snap_cleandata[0].df.columns
        
        if self.type == 'block':
            
            bl_number = []
            x = []
            y = []
            z = []
            
            
            # Compose long vectors of coordinates x,y,z
            
            for snap_bl in self.snap_cleandata:
                
                
                bl_number.append( snap_bl.num )
                
                x_bl = snap_bl.g.gx
                y_bl = snap_bl.g.gy
                z_bl = snap_bl.g.gz
                
                # Notice the order of output X,Y,Z !
                # https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html
                
                X,Y,Z = np.meshgrid( x_bl, y_bl, z_bl, indexing='ij' )
                
                # x,y,z are lists of nparrays
                
                x = x + np.ravel( X.T ).tolist()
                y = y + np.ravel( Y.T ).tolist()
                z = z + np.ravel( Z.T ).tolist()
                
            x = np.array( x ).ravel()
            y = np.array( y ).ravel()
            z = np.array( z ).ravel()
            
            GX = np.stack( [ x, y, z ] )
            
            GX_header = 'x y z '
            
            
            # Compose long vectors of solutions
                       
            sol_bl = np.zeros( (len(vars),self.n_cells_clean), dtype=np.float32 )
        
            pos_s = 0
            
            for snap_bl in self.snap_cleandata:
                
                pos_e = pos_s + snap_bl.npx*snap_bl.npy*snap_bl.npz
                
                sol_bl[:,pos_s:pos_e] = np.array( snap_bl.df.values ).T
                
                pos_s = pos_e
                
                
        elif self.type == 'slice':
            
            if self.slic_type == 'X':
                
                bl_number = []
                y = []
                z = []
                
                # Compose long vectors of coordinates y,z
                
                for snap_bl in self.snap_cleandata:
                    
                    bl_number.append(snap_bl.num)
                    
                    # GX => snap_bl[4] 
                    # For slice, only two coordinates vectors
                    
                    y_bl = snap_bl.g.gy
                    z_bl = snap_bl.g.gz
                    
                    Y, Z = np.meshgrid( y_bl, z_bl )
                    
                    y = y + np.ravel( Y ).tolist()
                    z = z + np.ravel( Z ).tolist()

                y = np.array( y ).ravel()
                z = np.array( z ).ravel()
                
                GX = np.stack( [ y, z ] )
                
                GX_header = 'y z '
                
                # Compose long vectors of solutions
                
                sol_bl = np.zeros( (len(vars),self.n_cells_clean), dtype=np.float32 )
                
                pos_s = 0 
                
                for snap_bl in self.snap_cleandata:
                                        
                    pos_e = pos_s + snap_bl.npx*snap_bl.npy*snap_bl.npz
                        
                    sol_bl[:,pos_s:pos_e] = np.array( snap_bl.df.values ).T
                        
                    pos_s = pos_e
                    
                    

            elif self.slic_type == 'W' or self.slic_type == 'Y':
                
                bl_number = []
                x = []
                z = []
                
                
                # Compose long vectors of coordinates x,z
                
                for snap_bl in self.snap_cleandata:
                    
                    
                    bl_number.append(snap_bl.num)
                    
                    x_bl = snap_bl.g.gx
                    z_bl = snap_bl.g.gz
                    
                    X, Z = np.meshgrid( x_bl, z_bl )
                    
                    x = x + np.ravel( X ).tolist()
                    z = z + np.ravel( Z ).tolist()

                x = np.array( x ).ravel()
                z = np.array( z ).ravel()
                
                GX = np.stack( [ x, z ] )
                
                GX_header = 'x z '
                

                # Compose long vectors of solutions
                
                sol_bl = np.zeros( (len(vars),self.n_cells_clean ), dtype=np.float32 )
                
                pos_s = 0 
                
                for snap_bl in self.snap_cleandata:
                                        
                    pos_e = pos_s + snap_bl.npx*snap_bl.npy*snap_bl.npz
                        
                    sol_bl[:,pos_s:pos_e] = np.array( snap_bl.df.values ).T
                        
                    pos_s = pos_e

            
            elif self.slic_type == 'Z':
                
                bl_number = []
                x = []
                y = []
                
                
                # Compose long vectors of coordinates x,y
                
                for snap_bl in self.snap_cleandata:
                    
                    
                    bl_number.append(snap_bl.num)
                    
                    x_bl = snap_bl.g.gx
                    y_bl = snap_bl.g.gy
                    
                    X, Y = np.meshgrid( x_bl, y_bl )
                    
                    x = x + np.ravel( X ).tolist()
                    y = y + np.ravel( Y ).tolist()

                x = np.array( x )
                y = np.array( y )
                
                GX = np.stack( [ x, y ] )
                
                GX_header = 'x y '
                

                # Compose long vectors of solutions
                
                sol_bl = np.zeros( (len(vars),self.n_cells_clean), dtype=np.float32 )
                
                pos_s = 0 
                
                for snap_bl in self.snap_cleandata:
                                        
                    pos_e = pos_s + snap_bl.npx*snap_bl.npy*snap_bl.npz
                        
                    sol_bl[:,pos_s:pos_e] = np.array( snap_bl.df.values ).T
                        
                    pos_s = pos_e
            
        
        # Return pandas dataframe
        
        
        df = pd.DataFrame( GX.T, columns=GX_header.strip().split() )
        
        df_sol = pd.DataFrame( sol_bl.T, columns=vars )
        
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
# >>> Get grid vectors                                           ( 7 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/05  - created
#
# Desc
#
#   - read a snapshot and return GX (DMD needs this)
#   - mainly because for 2D snapshots,we cannot have the grid info from grid.bin 
#
# ----------------------------------------------------------------------

    def get_grid_vectors( self, buff=3 ):
        
        """
        return GX (tall numpy array) for DMD
        """
        
        self.verbose = False
        
        self.read_snapshot()
        
        if self.type == 'block':
            
            x = []
            y = []
            z = []
        
            for snap_bl in self.snap_data:
                
                x_bl = snap_bl.g.gx[buff:-buff]
                y_bl = snap_bl.g.gy[buff:-buff]
                z_bl = snap_bl.g.gz[buff:-buff]
                
                # Notice the order of output X,Y,Z !
                # https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html

                X,Y,Z = np.meshgrid( x_bl, y_bl, z_bl, indexing='ij') 
                
                # x,y,z are lists of nparrays
                
                x.extend( np.ravel( X.T ).tolist() )
                y.extend( np.ravel( Y.T ).tolist() )
                z.extend( np.ravel( Z.T ).tolist() )
            
            # turn x into np.array, otherwise list cannot ravel
        
            GX = np.array([x,y,z]).T
        
            return GX
        
        
        elif self.type == 'slice':
            
            if self.slic_type == 'X':
                
                y=[]
                z=[]
                
                for snap_bl in self.snap_data:
                    
                    y_bl = snap_bl.g.gy[buff:-buff]
                    z_bl = snap_bl.g.gz[buff:-buff]
                    
                    Y,Z = np.meshgrid( y_bl, z_bl )
                    
                    y.extend( np.ravel( Y ).tolist() )
                    z.extend( np.ravel( Z ).tolist() )
                
                GX = np.array([y,z]).T
                
                return GX
            
            
            elif self.slic_type == 'W' or self.slic_type == 'Y':
                
                x=[]
                z=[]
                
                for snap_bl in self.snap_data:
                    
                    x_bl = snap_bl.g.gx[buff:-buff]
                    z_bl = snap_bl.g.gz[buff:-buff]
                    
                    X,Z = np.meshgrid( x_bl, z_bl )
                    
                    x.extend( np.ravel( X ).tolist() )
                    z.extend( np.ravel( Z ).tolist() )
                
                GX = np.array([x,z]).T
                
                return GX
            
            
            elif self.slic_type == 'Z':
                
                x=[]
                y=[]
                
                for snap_bl in self.snap_data:
                    
                    x_bl = snap_bl.g.gx[buff:-buff]
                    y_bl = snap_bl.g.gy[buff:-buff]
                    
                    X,Y = np.meshgrid( x_bl, y_bl )
                    
                    x.extend( np.ravel( X ).tolist() )
                    y.extend( np.ravel( Y ).tolist() )

                GX = np.array([x,y]).T
                
                return GX



# ----------------------------------------------------------------------
# >>> Get slice                                                   ( 8 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/07/31  - created
#
# Desc
#   
#   - get a 2D slice from a 3D snapshot
#   - inputs: slice normal to which axis; location of the slice
#
# ----------------------------------------------------------------------

    def get_slice( self, slic_type, loc:float, buff=3 ):
        
        """
        slic_type : 'X','Y' or 'Z'
        loc       : location of slice
        
        return    : 2D Snapshot object
        """
             
# ----- check if current snapshots 3D?
        
        if not self.type == 'block':
            raise TypeError("The snapshot to be sliced is not 3D")
        
        
# # ----- check if the grid file is available and read in grid
        
#         if not os.path.exists('inca_grid.bin'):
#             raise FileNotFoundError("Please check inca_grid.bin!")
        
#         else:
            
#             grid3d = GridData( 'inca_grid.bin' )
#             grid3d.verbose = True
#             grid3d.read_grid()

        if self.grid3d is None:
            raise ValueError("Please read in grid file first!")
        
        else: grid3d = self.grid3d   
        

# ----- select the blocks that intersect the plane
        # keep 1. intersected blocks' number  2. indexes of sliced location
        
        bl_intersect = []
        indx_slic    = []
        
        # loop over all blocks within this 3d snapshot
               
        for bl_data in self.snap_data:
            
            bl_num = bl_data.num      # block number starts from 1
            grd = grid3d.g[bl_num-1]  # index of grid starts from 0
            
            if slic_type == 'X':
                
                # find the bounding box of the block
                xmin = grd.lx0
                xmax = grd.lx1
                
                # if intersect, keep this block number and 
                # find the index of slice location
                if loc >= xmin and loc < xmax:
                    
                    bl_intersect.append( bl_num )
                    
                    istart = buff
                    iend   = grd.nx + buff
                    
                    for i in range(istart,iend):
                        
                        bnd1 = grd.gx[i] - 0.5*grd.hx[i]
                        bnd2 = grd.gx[i] + 0.5*grd.hx[i]
                        if loc >= bnd1 and loc < bnd2:
                            indx_slic.append(i)
                            break

            
            elif slic_type == 'Y':
                
                ymin = grd.ly0
                ymax = grd.ly1
               
                if loc >= ymin and loc < ymax:
                    bl_intersect.append( bl_num )

                    istart = buff
                    iend   = grd.ny + buff
                    
                    for i in range(istart,iend):
                        
                        bnd1 = grd.gy[i] - 0.5*grd.hy[i]
                        bnd2 = grd.gy[i] + 0.5*grd.hy[i]
                        if loc >= bnd1 and loc < bnd2:
                            indx_slic.append(i)
                            break
            
            
            elif slic_type == 'Z':
                
                zmin = grd.lz0
                zmax = grd.lz1
                
                if loc >= zmin and loc < zmax:
                    bl_intersect.append( bl_num )

                    istart = buff
                    iend   = grd.nz + buff
                    
                    for i in range(istart,iend):
                        
                        bnd1 = grd.gz[i] - 0.5*grd.hz[i]
                        bnd2 = grd.gz[i] + 0.5*grd.hz[i]
                        if loc >= bnd1 and loc < bnd2:
                            indx_slic.append(i)
                            break
        
        if len(bl_intersect) == 0:
            raise ValueError("No block is sliced! Check slice location.")
                    
        if self.verbose:
            print(f"sliced {len(bl_intersect)} blocks.\n")
            
            print("index of block  --  index of cut cell")
            for i in range(len(bl_intersect)):
                print(f"{bl_intersect[i]:15d} -- {indx_slic[i]:12d}")

# ----- init a slice snapshot 
        
        snap_2d = Snapshot()

        # Get sliced data from 3d snapshot and fill in 2d snapshot
        
        for bl_data in self.snap_data:
            
            bl_num = bl_data.num
            
            if bl_num in bl_intersect:
                
                grd = grid3d.g[bl_num-1]
                idx = indx_slic[bl_intersect.index(bl_num)]
                
                N1 = grd.nx + buff*2
                N2 = grd.ny + buff*2
                N3 = grd.nz + buff*2
                G1 = grd.gx
                G2 = grd.gy
                G3 = grd.gz
                
                vars   = bl_data.df.columns
                n_var  = len( vars )
                df_sol = deepcopy( bl_data.df.values )
                sol = np.array(df_sol).T.reshape(n_var,N3,N2,N1)
                
                if slic_type == 'X':                
                    
                    dims = [1,N2,N3]
                    GX = [G2,G3]
                    sol = sol[:,:,:,idx].reshape(n_var,N3*N2)
                                    
                if slic_type == 'Y':
                    
                    dims = [N1,1,N3]
                    GX = [G1,G3]
                    sol = sol[:,:,idx,:].reshape(n_var,N3*N1)
                    
                if slic_type == 'Z':
                    
                    dims = [N1,N2,1]
                    GX   = [G1,G2]
                    sol  = sol[:,idx,:,:].reshape(n_var,N2*N1)
                
                df_sol   = pd.DataFrame(sol.T,columns=vars)
                
                bl_slice = SnapBlock()
                bl_slice.fill_with_data( bl_num, dims, GX, df_sol, 'slice' )
                snap_2d.snap_data.append( bl_slice )            
        
        # fill in headers
        
        snap_2d.hformat      = self.hformat
        snap_2d.kind         = self.kind
        snap_2d.itstep       = self.itstep
        snap_2d.n_species    = self.n_species
        snap_2d.itime        = self.itime
        snap_2d.snap_lean    = self.snap_lean
        snap_2d.compressible = self.compressible
        snap_2d.snap_with_gx = self.snap_with_gx
        snap_2d.snap_with_tp = self.snap_with_tp
        snap_2d.snap_with_vp = self.snap_with_vp
        snap_2d.snap_with_cp = self.snap_with_cp
        snap_2d.snap_with_mu = self.snap_with_mu
        
        snap_2d.snap_with_cf = False
        snap_2d.snap_with_wd = self.snap_with_wd
        snap_2d.snap_with_bg = False  
        
        snap_2d.slic_type    = slic_type
        snap_2d.type         = 'slice'
        snap_2d.n_var        = self.n_var
        snap_2d.vars_name    = self.vars_name
        
        return snap_2d


# ----------------------------------------------------------------------
# >>> get slice dataframe (snapshot.bin                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def get_slice_df( self, slic_type, loc, buff=3 ):
        
        """
        slic_type : 'X','Y' or 'Z'
        loc       : location of slice
        
        return    : pandas dataframe of clean data(without ghost) on the slice
        """
        
        snapshot2d = self.get_slice( slic_type, loc, buff=buff )
        
        snapshot2d.drop_ghost( buff=3 )
        
        snapshot2d.assemble_block()
        
        return snapshot2d.df


# ----------------------------------------------------------------------
# >>> Get probed data                                       (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def get_probed_df( self, bl_list, G, indx_probed, probe_type, buff=3):
        
        """
        bl_list    : list of probed blocks.\n
        G          : corresponding GridData instance.\n
        indx_probed: list of indices tuple for each block.\n
        probe_type : 'X','Y','Z'\n
        
        return     : dataframe of probed results\n 
        """

# ----- extract probed line data on each bl.df

        indx_in_snap = [self.bl_nums.index(bl_num) for bl_num in bl_list]

        for i, bl_indx in enumerate( indx_in_snap ):
            
            bl_df = self.snap_data[bl_indx].df
            g     = G.g[bl_list[i]-1]
            
            indx  = indx_probed[i]
            
            npx   = g.nx + buff*2
            npy   = g.ny + buff*2
            npz   = g.nz + buff*2
            
            vars = bl_df.columns
            
            data_chunk = None
            
            for var in vars:
                
                data = np.array( bl_df[var] )
                data = data.reshape( npz, npy, npx )
                
                if   probe_type == 'X': data = data[indx[0],indx[1],:]
                elif probe_type == 'Y': data = data[indx[1],:,indx[0]]
                elif probe_type == 'Z': data = data[:,indx[1],indx[0]]
                
                if data_chunk is None: data_chunk = [data]
                else: data_chunk.append(data)
                
            data_chunk = np.array(data_chunk).T
            
            bl_df = pd.DataFrame(data_chunk,columns=vars)
            
# ------ match grids with data

            if probe_type == 'X': bl_df['x'] = g.gx
            if probe_type == 'Y': bl_df['y'] = g.gy
            if probe_type == 'Z': bl_df['z'] = g.gz
        
# ------ drop ghost cells

            self.snap_data[bl_indx].df = bl_df.iloc[buff:-buff]
            
# ------ concatenate all dataframe in all selected blocks

        df_probe = pd.concat( [self.snap_data[bl_indx].df 
                               for bl_indx in indx_in_snap] )
        
        df_probe.sort_values(by=[probe_type.lower()],inplace=True)
        
        df_probe.reset_index( drop=True, inplace=True )
        
        return df_probe


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/05/06  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def copy_var_from( self, snap_source, varnames:list):

        """
        snap_source: Snapshot instance where data is copied from
        varname    : name of variable to be copied
        """
        
        for snap_bl in self.snap_data:
            
            bl_num = snap_bl.num
            
            bl_src_indx = snap_source.bl_nums.index(bl_num)
            bl_src = snap_source.snap_data[bl_src_indx]
            
            for var in varnames:
                snap_bl.df[var] = bl_src.df[var]

        print(f"copied variables {varnames} from snapshot {snap_source.itstep} to snapshot {self.itstep}.\n")


# ----------------------------------------------------------------------
# >>> Write snapshot                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/02  - created
#
# Desc
#
#   - write a snapshot into binary file
#
# ----------------------------------------------------------------------

    def write_snapshot( self, filename ):
        
        """
        filename : filename of output snapshot
        """
        
        verbose = self.verbose
        
# ----- write snapshot header 
    
        # Default size
        sin = 4
        slg = 4
        sfl = 8
        
        with open( filename, 'wb' ) as f:
            
            pos = 0
            
            # header format
            
            write_int_bin( self.hformat, f, sin )
            
            # floating point precision format
            
            write_int_bin( self.kind, f, sin )
            
            # time step 
            
            write_int_bin( self.itstep, f, sin )
            
            # number of species
            
            write_int_bin( self.n_species, f, sin )
            
            pos += 4*sin

            # simulation time
            
            write_flt_bin( self.itime, f, sfl )
            
            pos += sfl
            
            # compose logicals
            
            if self.type == 'block':
                
                buf_log = [ self.snap_lean,
                            self.compressible,
                            self.snap_with_gx,
                            self.snap_with_tp,
                            self.snap_with_vp,
                            self.snap_with_cp,
                            self.snap_with_mu,
                            self.snap_with_wd ]
                
                write_log_bin( buf_log, f, slg )
                
                pos += slg*8
            
            elif self.type == 'slice':
                
                buf_log = [ self.snap_lean,
                            self.compressible,
                            self.snap_with_gx,
                            self.snap_with_tp,
                            self.snap_with_vp,
                            self.snap_with_cp,
                            self.snap_with_mu,
                            self.snap_with_cf,
                            self.snap_with_wd,
                            self.snap_with_bg ]

                write_log_bin( buf_log, f, slg )
                
                pos += slg*10
                                
            else: raise ValueError("snapshot type not supported.")
            
            if verbose: 
                print(f"Finish write snapshot header, pos = {pos}.")
            
# ----- write snapshot body

            for bl_data in self.snap_data:
                
                # block number
                
                write_int_bin( bl_data.num, f, sin )
                
                # dimensions of grids, N1, N2, N3
                
                write_int_bin( bl_data.npx, f, sin )
                write_int_bin( bl_data.npy, f, sin )
                write_int_bin( bl_data.npz, f, sin )
                
                # GX (G1,G2,G3 different dimension for different cases
                #  ,thus cannot be transformed into numpy array directly.)

                if self.type == 'block':
                    
                    write_flt_bin( bl_data.g.gx, f, self.kind )
                    write_flt_bin( bl_data.g.gy, f, self.kind )
                    write_flt_bin( bl_data.g.gz, f, self.kind )
                
                elif self.type == 'slice':
                    if bl_data.npx == 1:
                        write_flt_bin( bl_data.g.gy, f, self.kind )
                        write_flt_bin( bl_data.g.gz, f, self.kind )
                    if bl_data.npy == 1:
                        write_flt_bin( bl_data.g.gx, f, self.kind )
                        write_flt_bin( bl_data.g.gz, f, self.kind )
                    if bl_data.npz == 1:
                        write_flt_bin( bl_data.g.gx, f, self.kind )
                        write_flt_bin( bl_data.g.gy, f, self.kind )

                # sol(n_var, N3, N2, N1)
                
                write_flt_bin( np.array(bl_data.df.values).T, f, self.kind )
                

# ----------------------------------------------------------------------
# >>> Assign wall_dist field from wd_snap to a real snapshot      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/04/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def assign_wall_dist( self, wd_snap ):
        
        """
        wd_snap : Snapshot instance with wall distance field
        
        copy wall distance field from wd_snap to self.snap_data[].df
        """
        
        for snap_bl in self.snap_data:
            
            bl_num = snap_bl.num
            wd_df = wd_snap.snap_data[bl_num-1].df
            snap_bl.df['wd'] = wd_df['wd']


# ----------------------------------------------------------------------
# >>> compute bubble volume                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/04/22  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_bubble_volume( self, G:GridData, cc_df=None, 
                               roughwall=False, buff=3 ):
        
        """
        G     : GridData instance
        cc_df : cutcell dataframe from cutcells_setup.dat
        
        return: separation bubble volume
        Need data chunk with u ready.
        G should contain cell volume.
        wd (wall distance) should be contained in self.snap_data[num-1].df
        """
        
        vol_bubble = 0.0
        
        for snap_bl in self.snap_data:
            
            bl_num = snap_bl.num
            g = G.g[bl_num-1]
            
            npx = g.nx + buff*2
            npy = g.ny + buff*2
            npz = g.nz + buff*2
            
            data_df = snap_bl.df
            
            u = np.array(data_df['u']).reshape(npz,npy,npx)
            
            identifier = u < 0.0
            identifier = identifier*1.0   # convert to float
            
            vol = g.vol
            
            if roughwall:
            
                temp_df = cc_df[cc_df['block_number'] == bl_num]
                wall_dist = np.array( snap_bl.df['wd'] )
                g.assign_vol_fra( df=temp_df, wall_dist=wall_dist )
            
            else:
                
                g.assign_vol_fra()
                
#            print(np.shape(vol),np.shape(identifier),np.shape(g.vol_fra.T))
            
            vol_bubble_block = vol*identifier*(g.vol_fra.T)
                
            vol_bubble += np.sum(vol_bubble_block[buff:-buff,buff:-buff,buff:-buff])
        
        self.vol_bubble = vol_bubble
        
        return vol_bubble    


# ----------------------------------------------------------------------
# >>> compute separatio bubble volume from PDF                    (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/05/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_bubble_volume_pdf( self, G:GridData, cc_df=None,
                                   roughwall=False, opt=1, buff=3 ):
        
        """
        compute separation bubble volume from PDF of separation.
        
        Need the pdf_sep ready in self.snap_data[].df
        """
        
        vol_bubble = 0.0
        
        for snap_bl in self.snap_data:
            
            bl_num = snap_bl.num
            g = G.g[bl_num-1]
            
            npx = g.nx + buff*2
            npy = g.ny + buff*2
            npz = g.nz + buff*2
            
            data_df = snap_bl.df
            
            pdf_sep = np.array(data_df['pdf_sep']).reshape(npz,npy,npx)
            vol = g.vol
            
            if opt == 1:
                
                identifier = pdf_sep > 0.5
                identifier = identifier*1.0   # convert to float
                
                if roughwall:
                    
                    temp_df = cc_df[cc_df['block_number'] == bl_num]
                    wall_dist = np.array( snap_bl.df['wd'] )
                    g.assign_vol_fra( df=temp_df, wall_dist=wall_dist )
                
                else: g.assign_vol_fra()
            
                vol_bubble_block = vol*identifier*(g.vol_fra.T)
                vol_bubble += np.sum(vol_bubble_block[buff:-buff,buff:-buff,buff:-buff])
            
            elif opt == 2:
                
                if roughwall:
                        
                    temp_df = cc_df[cc_df['block_number'] == bl_num]
                    wall_dist = np.array( snap_bl.df['wd'] )
                    g.assign_vol_fra( df=temp_df, wall_dist=wall_dist )    
        
                else: g.assign_vol_fra()
                
                vol_bubble_block = vol*pdf_sep*(g.vol_fra.T)
                vol_bubble += np.sum(vol_bubble_block[buff:-buff,buff:-buff,buff:-buff])
        
        self.vol_bubble = vol_bubble
        
        return vol_bubble
    

# ----------------------------------------------------------------------
# >>> Write snapshot into tecplot szplt format                   (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/29  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def write_szplt( self, filename, vars=None, block_list=None, buff=3 ):
        
        """
        filename : filename of output snapshot
        
        tbd: output wall distance as well
        """

   
# ----- check block list, if None, write all blocks

        if block_list is None:
            block_list = self.bl_nums
        
# ----- drop ghost cells

        self.drop_ghost( block_list=block_list, buff=buff )        
        
# ----- check the variables to be written

        if vars is None:
            vars = self.snap_cleandata[0].df.columns.tolist()
            
# ----- setupt tecplot file
        
        tp.new_layout()
        frame = tp.active_frame()
        dataset = frame.create_dataset('snapshot',['x','y','z'] + vars)
        
        for i in range( len(self.snap_cleandata) ):
            
            bl_data = self.snap_cleandata[i]
            
            npx = int(bl_data.npx)
            npy = int(bl_data.npy)
            npz = int(bl_data.npz)
            
            zone = dataset.add_ordered_zone( f'bl_{bl_data.num:05d}',
                                             (npz, npy, npx),
                                             dtypes=FieldDataType.Float )

            xx,yy,zz = np.meshgrid( bl_data.g.gx, bl_data.g.gy, bl_data.g.gz, indexing='ij' )

            # tecplot needs x,y,z in C order
            
            zone.values('x')[:] = xx.ravel()
            zone.values('y')[:] = yy.ravel()
            zone.values('z')[:] = zz.ravel()
            
            for var in vars:
                zone.values(var)[:] = np.array(bl_data.df[var]).reshape(npz,npy,npx).T.ravel()
            
        tp.data.save_tecplot_szl( filename )
        print(f"Finish writing snapshot to {filename}.")
        

# ----------------------------------------------------------------------
# >>> create vtk multiblock dataset (snapshot)                   (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def create_vtk_multiblock( self, vars=None, block_list=None, buff=3, mode='symmetry' ):
        
        """
        write snapshot into vtm file (multiblock vtk)\n        
        vars       : list of variables to be written\n
        block_list : block numbers of which blocks to be written\n
        buff       : number of ghost layers\n
        mode       : 'symmetry' or 'oneside'
        """

# ----- check block list, if None, write all blocks

        if block_list is None:
            block_list = self.bl_nums
            
# ----- check if grid data is ready
        
        if self.grid3d is None:
            raise ValueError("Please read in grid data first!")
        
        else: 
            G = self.grid3d
        
# ----- drop ghost cells

        self.drop_ghost( block_list=block_list, buff=buff, mode=mode )
        
        if mode == 'symmetry':  buffl = buff; buffr = buff
        elif mode == 'oneside': buffl = buff; buffr = buff-1

# ----- check the variables to be written

        if vars is None:
            vars = self.snap_cleandata[0].df.columns.tolist()
            
# ----- setup vtk file
        
        vtk_blocks = list()
        
        for snap_bl in self.snap_cleandata:
            
            if snap_bl.num not in block_list:
                continue
            
            bl_num = snap_bl.num
            g = G.g[bl_num-1]
            
            px = g.px[buffl:-buffr]
            py = g.py[buffl:-buffr]
            pz = g.pz[buffl:-buffr]
            
            # modify the grid points arrays based on snapshot type
            if self.type == 'block':
                pass
            elif self.type == 'slice':
                if self.slic_type   == 'X': px = np.array([0.0])
                elif self.slic_type == 'Z': pz = np.array([0.0])
                elif self.slic_type == 'Y' or self.slic_type == 'W': py = np.array([0.0])
            
            # build one vtk block
            bl_vtk = create_3d_vtkRectilinearGrid( px, py, pz )
            
            for var in vars:
                
                var_data = np.array(snap_bl.df[var])
                
                if len(var_data) != snap_bl.npx*snap_bl.npy*snap_bl.npz:
                    raise ValueError(f"Data length not match for variable {var} in block {bl_num}.")
                elif len(var_data) == 0:
                    raise ValueError(f"Data length is zero for variable {var} in block {bl_num}.")
                
                bl_vtk = add_var_vtkRectilinearGrid( bl_vtk, var, var_data )
                
            vtk_blocks.append( bl_vtk )
        
        # build the multiple blocks dataset
        dataset = create_multiblock_dataset(vtk_blocks)

        return dataset


# ----------------------------------------------------------------------
# >>> write snapshot into vtm (multiblock vtk) file            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/06  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def write_vtm( self, filename, vars, block_list=None, buff=3 ):
        
        """
        write snapshot into vtm file (multiblock vtk)
        
        filename : filename of output snapshot
        vars     : list of variables to be written
        """

# ----- build the multiple blocks dataset
        dataset = self.create_vtk_multiblock(vars, block_list=block_list, buff=buff)
        
        write_vtm_file( filename, dataset )
        

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

    # test writing snapshot into vtm
    
    
    test_dir  = '/media/wencan/Expansion/temp/smooth_awallrs_test/snapshots/snapshot_00773347'
    snapshot_file = test_dir + '/snapshot.bin'
    grid_file = '/media/wencan/Expansion/temp/smooth_awallrs_test/results/inca_grid.bin'
    vars = ['u','v','w','p','T']
    
    os.chdir( test_dir )
    
    G = GridData( grid_file )

    G.read_grid()
    
    snapshot1 = Snapshot( snapshot_file )
    snapshot1.grid3d = G    
    snapshot1.verbose = False
    
    with timer('read one snapshot '):
        
        snapshot1.read_snapshot( var_read = vars )
        
        print(snapshot1.snap_data[0].df)
        print(snapshot1.snap_data[1].g.gx)
        print(snapshot1.snap_data[1].g.gy)

    with timer('write into vtm'):
        
        snapshot1.write_vtm( 'test.vtm', vars )

    # check y plane slice
    
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
          
    pass 

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
