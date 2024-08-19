# -*- coding: utf-8 -*-
'''
@File    :   grid.py
@Time    :   2023/02/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   GridData is a class referring to a grid data binary file.
             GridBlock is the data structure class within the range of a block.
'''


import os
import numpy             as np
import pandas            as pd

from   .io_binary        import read_int_bin
from   .io_binary        import read_flt_bin
from   .io_binary        import read_log_bin
from   .io_binary        import read_char_bin
from   .tools            import is_above_wavywall
from   .tools            import if_overlap_3d
from   .tools            import mean_of_list
from   .tools            import point_in_box

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
    
    def __init__(self, grid_file=None ):
        
        # directory to the grid binary file
        self.grid_file = grid_file
        
        # grid file size
        self.fsize = 0
        
        # file pointer position
        self.pos = 0
        
        # total number of grids with 3 layers ghost cells 
        self.n_total = 0
        
        # number of grids without ghost cells
        self.n_real  = 0
        
        # GridBlock List
        self.g = list()
        
        # verbose ?
        self.verbose = False



# ----------------------------------------------------------------------
# >>> Read Grid                                                (Nr.)
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
#   - read the grid file 
#
# ----------------------------------------------------------------------

    def read_grid( self ):

        # grid file size
        self.fsize = os.stat( self.grid_file ).st_size
        
        with open( self.grid_file, 'rb' ) as f:
            
            self.read_grid_header( f )
            
            self.read_grid_body( f )
        
        print( f"finish read grid file ...{self.grid_file[-50:]} \n" )



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
        
        self.n_bl = read_int_bin( file.read(sin), sin )
        self.pos += sin
        
        self.pos += 4
        
        self.header_size = self.pos
        

        if self.verbose:
            
            print(f"grid_file_format: {self.grid_file_format}")
            print(f"grid with solver: {self.grid_with_solver}")
            print(f"number of blocks: {self.n_bl}")

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
        # number start from 1
        
        i = 1
        
#        self.g.append( GridBlock( file, 0,self.grid_with_solver))
        
        while not end_of_file:
            
            # read in block grid one by one
            # input: file, index, if grid_with_solver
                        
            self.g.append(GridBlock(file,i,self.grid_with_solver,self.verbose))
            
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
#    - looks like useless, no where use it.
#    - low Re number case, there are 8 blocks in spanwise direction, 
#      this method helps to find the block group(sharing same lx0,ly0)
#
# ----------------------------------------------------------------------

    def get_sorted_groups( self ):
    
        """
        based on lx0,ly0,lz0, sorting block grids index (not Blockgrid!)
        same x-y location blocks form a block group
        """        
        if len(self.g) == 0:
            raise ValueError('Please read grid blocks first!')
        
        lx0 = list()
        ly0 = list()
        lx1 = list()
        ly1 = list()
        lz0 = list()
        bl_index = list()
        
        # get lx0,ly0,lz0 and bl_index lists
        for i in range(self.n_bl):
            
            lx0.append(self.g[i].lx0)
            ly0.append(self.g[i].ly0)
            lx1.append(self.g[i].lx1)
            ly1.append(self.g[i].ly1)
            lz0.append(self.g[i].lz0)
            bl_index.append(self.g[i].num)
        
        df = pd.DataFrame(lx0,columns=['lx0'])
        df['ly0'] = ly0
        df['lx1'] = lx1
        df['ly1'] = ly1
        df['lz0'] = lz0
        df['bl_index'] = bl_index
        
        # sort based on lx0,ly0,lz0
        df_sorted = df.sort_values(by=['lz0','ly0','lx0'])
        
        # group pd rows and aggregate lz0 and bl_index into lists
        self.grouped  = df_sorted.groupby(by=['ly0','lx0']).agg(list).\
                                  reset_index()

        self.grouped['lx1'] = self.grouped['lx1'].apply(mean_of_list)
        self.grouped['ly1'] = self.grouped['ly1'].apply(mean_of_list)


# ----------------------------------------------------------------------
# >>> Get Selected Block Grids                                 ( 1-4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/06  - created
#
# Desc
#
# - input a bbox[xmin,xmax,ymin,ymax,zmin,zmax]
# - record which block grids groups are partly or fully overlapped
#   with rect1 (mode=='overlap); or within bounding box(mode='within')
#
# ----------------------------------------------------------------------

    def select_blockgrids( self, bbox, mode='overlap' ):
        
        """
        bbox : [xmin,xmax,ymin,ymax,zmin,zmax]   \n
        mode : 'overlap' or 'within'
        """
        
        selected_bls = []
        
        for i in range( len(self.g) ):
            
            bl_num = self.g[i].num
            
            lx0 = self.g[i].lx0
            lx1 = self.g[i].lx1
            ly0 = self.g[i].ly0
            ly1 = self.g[i].ly1
            lz0 = self.g[i].lz0
            lz1 = self.g[i].lz1
            
            bbox2 = [ lx0, lx1, ly0, ly1, lz0, lz1 ]
            
            if mode == 'overlap':
                
                if if_overlap_3d( bbox, bbox2 ):
                    
                    selected_bls.append( bl_num )
                    
            elif mode == 'within':
                
                if (    lx0 >= bbox[0] and lx1 <= bbox[1]
                    and ly0 >= bbox[2] and ly1 <= bbox[3]
                    and lz0 >= bbox[4] and lz1 <= bbox[5] ):
                    
                    selected_bls.append( bl_num )
            
        return selected_bls



# ----------------------------------------------------------------------
# >>> select sliced block grids                                 (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def select_sliced_blockgrids( self, slic_type, loc, bbox=None, buff=3 ):
        
        """
        input:\n
        slic_type : 'X', 'Y', or 'Z' ( normal direction of slice) \n
        loc       : coordinate \n
        bbox      : [xmin,xmax,ymin,ymax,zmin,zmax] \n
        
        return : list of selected bl_nums, list of slice indexes on each bl
        """
        
        selected_bls = []
        indx_slic    = []
        
        if bbox is None: 
            min = float('-inf')
            max = float('inf')
            bbox = [min,max,min,max,min,max]
        
        within_bbox_bls = self.select_blockgrids( bbox, mode='within' )
        
        
        for grd in self.g:
            
            bl_num = grd.num
            
            lx0 = grd.lx0
            lx1 = grd.lx1
            ly0 = grd.ly0
            ly1 = grd.ly1
            lz0 = grd.lz0
            lz1 = grd.lz1
            
            if slic_type == 'X':
                
                if (lx0 <= loc < lx1) and (bl_num in within_bbox_bls):
                    
                    selected_bls.append( bl_num )
                    istart = buff
                    iend   = grd.nx + buff
                    
                    for i in range( istart, iend ):
                        
                        bnd1 = grd.gx[i] - 0.5*grd.hx[i]
                        bnd2 = grd.gx[i] + 0.5*grd.hx[i]
                        if (bnd1 <= loc < bnd2):
                            indx_slic.append(i)
                            break
            
            elif slic_type == 'Y':
                
                if (ly0 <= loc < ly1) and (bl_num in within_bbox_bls):
                    
                    selected_bls.append( bl_num )
                    istart = buff
                    iend   = grd.ny + buff
                    
                    for i in range( istart, iend ):
                        
                        bnd1 = grd.gy[i] - 0.5*grd.hy[i]
                        bnd2 = grd.gy[i] + 0.5*grd.hy[i]
                        if (bnd1 <= loc < bnd2):
                            indx_slic.append(i)
                            break                    
                    
            elif slic_type == 'Z':
                
                if (lz0 <= loc < lz1) and (bl_num in within_bbox_bls):
                    
                    selected_bls.append( bl_num )

                    istart = buff
                    iend   = grd.nz + buff
                    
                    for i in range( istart, iend ):
                        
                        bnd1 = grd.gz[i] - 0.5*grd.hz[i]
                        bnd2 = grd.gz[i] + 0.5*grd.hz[i]
                        if (bnd1 <= loc < bnd2):
                            indx_slic.append(i)
                            break
         
        if len(indx_slic) != len(selected_bls):
            raise ValueError("Number of slices index != number of blocks")


        return selected_bls, indx_slic


# ----------------------------------------------------------------------
# >>> Select probed block grids                                  (Nr.)
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

    def select_probed_blockgrids( self, probe_type, loc, bbox=None, buff=3):
        
        """
        probe_type: 'X', 'Y', or 'Z' (parallel direction of probing) \n
        loc       : [z,y],[x,z], or [x,y] (follow conventional order)\n
        bbox      : [xmin,xmax,ymin,ymax,zmin,zmax] \n
        
        return: a list of blocks and a list of indices tuple
        """
        
        selected_bls = []
        indx_probe   = []

        if bbox is None: 
            min = float('-inf')
            max = float('inf')
            bbox = [min,max,min,max,min,max]
        
        within_bbox_bls = self.select_blockgrids( bbox, mode='within' )
        
        
        for grd in self.g:
            
            bl_num = grd.num
            lx0 = grd.lx0
            lx1 = grd.lx1
            ly0 = grd.ly0
            ly1 = grd.ly1
            lz0 = grd.lz0
            lz1 = grd.lz1
            
            if probe_type == 'X':
                
                if ((lz0 <= loc[0] < lz1) and (ly0 <= loc[1] < ly1) and
                    (bl_num in within_bbox_bls)):
                    
                    selected_bls.append( bl_num )
                    
                    for i in range( buff, grd.nz + buff ):
                        
                        bnd1 = grd.gz[i] - 0.5*grd.hz[i]
                        bnd2 = grd.gz[i] + 0.5*grd.hz[i]
                        if (bnd1 <= loc[0] < bnd2):
                            break

                    for j in range( buff, grd.ny + buff ):
                        
                        bnd1 = grd.gy[j] - 0.5*grd.hy[j]
                        bnd2 = grd.gy[j] + 0.5*grd.hy[j]
                        if (bnd1 <= loc[1] < bnd2):
                            break
                    
                    indx_probe.append((i,j))
            
            
            elif probe_type == 'Y':
                
                if ((lx0 <= loc[0] < lx1) and (lz0 <= loc[1] < lz1) and
                    (bl_num in within_bbox_bls)):
                    
                    selected_bls.append( bl_num )
                    
                    for i in range( buff, grd.nx + buff ):
                        
                        bnd1 = grd.gx[i] - 0.5*grd.hx[i]
                        bnd2 = grd.gx[i] + 0.5*grd.hx[i]
                        if (bnd1 <= loc[0] < bnd2):
                            break

                    for j in range( buff, grd.nz + buff ):
                        
                        bnd1 = grd.gz[j] - 0.5*grd.hz[j]
                        bnd2 = grd.gz[j] + 0.5*grd.hz[j]
                        if (bnd1 <= loc[1] < bnd2):
                            break
                    
                    indx_probe.append((i,j)) 
                                  
                    
            elif probe_type == 'Z':
                
                if ((lx0 <= loc[0] < lx1) and (ly0 <= loc[1] < ly1) and 
                    (bl_num in within_bbox_bls)):
                    
                    selected_bls.append( bl_num )
                    
                    for i in range( buff, grd.nx + buff ):
                        
                        bnd1 = grd.gx[i] - 0.5*grd.hx[i]
                        bnd2 = grd.gx[i] + 0.5*grd.hx[i]
                        if (bnd1 <= loc[0] < bnd2):
                            break

                    for j in range( buff, grd.ny + buff ):
                        
                        bnd1 = grd.gy[j] - 0.5*grd.hy[j]
                        bnd2 = grd.gy[j] + 0.5*grd.hy[j]
                        if (bnd1 <= loc[1] < bnd2):
                            break
                    
                    indx_probe.append([i,j])
         
        if len(indx_probe) != len(selected_bls):
            raise ValueError("Number of slices index != number of blocks")

        return selected_bls, indx_probe
            

# ----------------------------------------------------------------------
# >>> Find probe index                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/08  - created
#
# Desc
#     - maybe used in finding solution's index
# ----------------------------------------------------------------------

    def find_probe_index( self, xyz, buff=3 ):
        
        """
        xyz: [x,y,z] of the given point
        return: bl_num, [i,j,k]
        Note: i,j,k are the index when buffer layers are included
        """

        for grd in self.g:
            
            if (grd.lx0 <= xyz[0] < grd.lx1 and
                grd.ly0 <= xyz[1] < grd.ly1 and
                grd.lz0 <= xyz[2] <= grd.lz1):
                
                i,j,k = grd.point_index( xyz, buff )
                bl_num = grd.num
                
                break
            
            
                
        return bl_num, [i, j, k]
                

# ----------------------------------------------------------------------
# >>> find probe xyz                                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/17  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def find_probe_xyz( self, xyz, buff=3 ):
        
        """"
        xyz: [x,y,z] of the given point
        """
        bl_num, indx = self.find_probe_index( xyz, buff=buff )
        
        g = self.g[bl_num-1]
        x = g.gx[indx[0]]
        y = g.gy[indx[1]]
        z = g.gz[indx[2]]
        
        return [x,y,z]


# ----------------------------------------------------------------------
# >>> compute cell volume                                  (Nr.)
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

    def cell_volume( self ):
        
        """
        compute the volume of each cell
        """
        
        for gblock in self.g:
            
            gblock.cell_volume()
        

# ----------------------------------------------------------------------
# >>> compute cell point (coordinates of vertice)                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_point( self ):
        
        for gblock in self.g:
            gblock.compute_point()


# ----------------------------------------------------------------------
# >>> group blocks by their range           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def group_by_range( self, grp_type, block_list=None ):
        
        """
        group blocks by their range
        
        grp_type : 'xy'('yx'), 'yz'('zy'), or 'zx'('xz)
        block_list : list of block numbers that to be grouped.
        
        return : a list of grouped blocks numbers.
        """

# ----- if block_list is None, then group all blocks        
        if block_list is None:
            block_list = [i for i in range(1, self.n_bl+1)]

# ----- read in grid blocks into pandas dataframe

        df = pd.DataFrame(columns=['lx0','lx1','ly0','ly1','lz0','lz1','bl_num'])
        
        for bl_num in block_list:
            
            bl = self.g[bl_num-1]
            df.loc[len(df)] = [bl.lx0, bl.lx1, bl.ly0, bl.ly1, bl.lz0, bl.lz1, bl_num]
        
        df['bl_num'] = df['bl_num'].astype(int)
            
# ----- group blocks by their range

        if grp_type == 'xy' or grp_type == 'yx':
            grouped = df.groupby(['lx0','lx1','ly0','ly1']).agg(list).reset_index()
            grouped.sort_values(by=['ly0','lx0'], inplace=True)
            
        elif grp_type == 'yz' or grp_type == 'zy':
            grouped = df.groupby(['ly0','ly1','lz0','lz1']).agg(list).reset_index()
            grouped.sort_values(by=['lz0','ly0'], inplace=True)
            
        elif grp_type == 'zx' or grp_type == 'xz':
            grouped = df.groupby(['lz0','lz1','lx0','lx1']).agg(list).reset_index()
            grouped.sort_values(by=['lx0','lz0'], inplace=True)
        
        grouped_block_list = grouped['bl_num'].tolist()
        
        # pd.set_option('display.max_rows', None)
        # print(grouped)
        
        return grouped_block_list
    

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

class GridBlock:
        
    def __init__( self, file=None , num=None , grid_with_solver=None, verbose=False ):
        
        """
        a void GridBlock can be initialized by given nothing or ...
        
        file: grid file
        num : index of GridBlock
        grid_with_solver : if grid with solver
        verbose : if print out info
        """
        
        if file is None:
            self.gx = float()
            self.dx = float()
            self.hx = float()
            self.wx = float()
            
            self.gy = float()
            self.dy = float()
            self.hy = float()
            self.wy = float()

            self.gz = float()
            self.dz = float()
            self.hz = float()
            self.wz = float()
            
            self.vol_fra = None
            self.vol = None
        
        else: 
            self._init_from_file( file, num, grid_with_solver, verbose )
            self.vol_fra = None
            self.vol = None

    def _init_from_file( self, file, num, grid_with_solver, verbose ):
        
        sin = 4
        sfl = 8
        slg = 4
        
        self.verbose = verbose
        self.grid_with_solver = grid_with_solver
        
        len_specname = 19      # should be 19, after INCA 5b66ccd, 7th Feb 2023 
                               # 15 before INCA 5b66ccd, 7th Feb 2023
        
        # index of GridBlock, like first block's grid, second ...
        self.num = num
        if self.verbose: print( 'grid number is', self.num )
        
        # size of this GridBlock (in bytes)
        self.size = 0
        
        # pass spaces before 'write'
        file.read(4)
        self.size += 4
        
        # read grid dimension 
        # nx, ny, nz do not include buffer cells, 
        # instead, npx,npy,npz should include buffer cells 
        # self.np number of cells in a block
        self.nx = read_int_bin( file.read(sin), sin )
        self.ny = read_int_bin( file.read(sin), sin )
        self.nz = read_int_bin( file.read(sin), sin )
        self.np = (self.nx+6)*(self.ny+6)*(self.nz+6)
        
        self.size += 3*sin
        if self.verbose: 
            print( 'nx = ', self.nx )
            print( 'ny = ', self.ny )
            print( 'nz = ', self.nz )
        
        # read grid bounding box corner location ( no ghost cells )
        self.lx0 = read_flt_bin( file.read(sfl), sfl )
        self.lx1 = read_flt_bin( file.read(sfl), sfl )
        self.ly0 = read_flt_bin( file.read(sfl), sfl )
        self.ly1 = read_flt_bin( file.read(sfl), sfl )
        self.lz0 = read_flt_bin( file.read(sfl), sfl )
        self.lz1 = read_flt_bin( file.read(sfl), sfl )
        
        self.size += 6*sfl
        if self.verbose:
            print( f'lx0 =  {self.lx0:10.5f} lx1 = {self.lx1:10.5f}' )
            print( f'ly0 =  {self.ly0:10.5f} ly1 = {self.ly1:10.5f}' )
            print( f'lz0 =  {self.lz0:10.5f} lz1 = {self.lz1:10.5f}' )
        
        # read parameters of grid: 3 components in each direction
        self.paramx = read_flt_bin( file.read(3*sfl), sfl )
        self.paramy = read_flt_bin( file.read(3*sfl), sfl )
        self.paramz = read_flt_bin( file.read(3*sfl), sfl )
        
        self.size += 9*sfl
        if self.verbose: print( 'paramx = ', self.paramx )
        
        # read shapes of grids: char(len=10) for each direction
        self.shapex = read_char_bin( file.read(10) )
        self.shapey = read_char_bin( file.read(10) )
        self.shapez = read_char_bin( file.read(10) )
        
        self.size += 30
        if self.verbose: print( 'shapex = ', self.shapex )
        
        # read boundary values: char(len=10)*2 for each direction
        self.bx1 = read_char_bin( file.read(10) )
        self.bx2 = read_char_bin( file.read(10) )
        self.by1 = read_char_bin( file.read(10) )
        self.by2 = read_char_bin( file.read(10) )
        self.bz1 = read_char_bin( file.read(10) )
        self.bz2 = read_char_bin( file.read(10) )
        
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
        # after INCA 5b66ccd, 7th Feb 2023, fluid names are 19 characters
        self.fluidx1 = read_char_bin( file.read(len_specname) )
        self.fluidx2 = read_char_bin( file.read(len_specname) )
        self.fluidy1 = read_char_bin( file.read(len_specname) )
        self.fluidy2 = read_char_bin( file.read(len_specname) )
        self.fluidz1 = read_char_bin( file.read(len_specname) )
        self.fluidz2 = read_char_bin( file.read(len_specname) )
        
        self.size += len_specname*6
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
        
        self.gx = read_flt_bin( file.read( (self.nx+6)*sfl ), sfl )
        self.dx = read_flt_bin( file.read( (self.nx+6)*sfl ), sfl )
        self.hx = read_flt_bin( file.read( (self.nx+6)*sfl ), sfl )
        self.wx = read_flt_bin( file.read( (self.nx+6)*sfl ), sfl )
        
        self.gy = read_flt_bin( file.read( (self.ny+6)*sfl ), sfl )
        self.dy = read_flt_bin( file.read( (self.ny+6)*sfl ), sfl )
        self.hy = read_flt_bin( file.read( (self.ny+6)*sfl ), sfl )
        self.wy = read_flt_bin( file.read( (self.ny+6)*sfl ), sfl )

        self.gz = read_flt_bin( file.read( (self.nz+6)*sfl ), sfl )
        self.dz = read_flt_bin( file.read( (self.nz+6)*sfl ), sfl )
        self.hz = read_flt_bin( file.read( (self.nz+6)*sfl ), sfl )
        self.wz = read_flt_bin( file.read( (self.nz+6)*sfl ), sfl )
        
        self.size += ( self.nx+self.ny+self.nz + 18 ) * sfl * 4
        if self.verbose: 
            print( 'self.gx = ', self.gx )
            print( 'self.gy = ', self.gy )
            print( 'self.gz = ', self.gz )
        
        # read extra spaces after 'write'
        file.read(4)
        self.size += 4



# ----------------------------------------------------------------------
# >>> compute corner point                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_point( self, buff=3 ):
        
        """
        compute corner points location
        
        return self.px,self.py,self.pz with length of nx(ny,nz) + buff*2 + 1
        """
        
        self.px = np.zeros( self.nx + buff*2 + 1 )
        self.py = np.zeros( self.ny + buff*2 + 1 )
        self.pz = np.zeros( self.nz + buff*2 + 1 )
        
        self.px[0] = self.gx[0] - 0.5*self.hx[0]
        self.py[0] = self.gy[0] - 0.5*self.hy[0]
        self.pz[0] = self.gz[0] - 0.5*self.hz[0]
        
        self.px[1:] = self.gx + 0.5*self.hx
        self.py[1:] = self.gy + 0.5*self.hy
        self.pz[1:] = self.gz + 0.5*self.hz
        

# ----------------------------------------------------------------------
# >>> Define Cut Cell                                          ( 2-1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/05  - created
# 2023/10/02  - revised, using wall distance field
#
# Desc
#
# - using dataframe containing i,j,k and volume fraction
# - 1. assign cut-cells volume fractions
# - 2. for cells, which is above wall and is not cut cell, vol_fra = 1.0
# ----------------------------------------------------------------------

    def assign_vol_fra( self, df=None, wall_dist=None ):
        
        """
        df        : cutcell info dataframe\n
        wall_dist : array of wall distance of a block, 1D, order as Snapshot\n
        
        return: self.vol_fra
        """
        
        if (df is not None) and len(df) > 0:
                
            # check block number
            if df['block_number'].iloc[0] != self.num :
                raise ValueError("block number does not match!")
            
            # take out data from df
            i = np.array( df['i'] )
            j = np.array( df['j'] )
            k = np.array( df['k'] )
            vol = np.array( df['vol'] )
            
            wall_dist = wall_dist.reshape(self.nz+6,self.ny+6,self.nx+6)
            wall_dist = wall_dist.T
            
            self.vol_fra = np.zeros( shape=(self.nx+6,self.ny+6,self.nz+6) )
            
            # for cells that above wall and not cut cell, set vol_fra = 1
            self.vol_fra[np.where(wall_dist >= 0.0)] = 1.0
            
            # df['i'], df['j'], df['k'] all count from 1 like Fortran
            # well in python, array count from 0.
            for index, value in enumerate( vol ):
                
                self.vol_fra[i[index]-1,j[index]-1,k[index]-1] = value
            
    #        print(f"Block {self.num} cut cell volume fractions are assigned.")

        # df is None: smooth wall ; len(df) ==0 : no cut cell
        elif (df is None) or len(df) == 0:
            self.vol_fra = np.ones( shape=(self.nx+6,self.ny+6,self.nz+6) )


# ----------------------------------------------------------------------
# >>> find point index                                   (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def point_index( self, xyz, buff=3 ):
        
        """
        xyz: list of x,y,z coordinates of the given point
        
        return: i,j,k
        """

        # check if the point is within the block
        
        bbox = [self.lx0, self.ly0, self.lz0, self.lx1, self.ly1, self.lz1]
        
        if not point_in_box( xyz, bbox ):
            raise ValueError(f"The point{xyz} is not within the block {self.num}")
                
        else:

            for i in range( buff, self.nx + buff ):
                if (self.gx[i]-0.5*self.hx[i] <= xyz[0] < self.gx[i]+0.5*self.hx[i]):
                    break
            
            for j in range( buff, self.ny + buff ):
                if (self.gy[j]-0.5*self.hy[j] <= xyz[1] < self.gy[j]+0.5*self.hy[j]):
                    break
            
            for k in range( buff, self.nz + buff ):
                if (self.gz[k]-0.5*self.hz[k] <= xyz[2] < self.gz[k]+0.5*self.hz[k]):
                    break

        return i,j,k


# ----------------------------------------------------------------------
# >>> compute cell volume                                         (Nr.)
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

    def cell_volume( self, buff=3 ):
        
        """
        compute cell volume, store results in self.vol.
        with buffer cells included
        """
        
        hx = self.hx; npx = self.nx + 2*buff
        hy = self.hy; npy = self.ny + 2*buff
        hz = self.hz; npz = self.nz + 2*buff
        
        self.vol = np.outer( hz, np.outer( hy, hx) ).reshape(npz,npy,npx)

# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    filename = '/media/wencanwu/Seagate Expansion Drive1/temp/231124/results/inca_grid.bin'
    
    grd = GridData( filename )
    
#    grd.grid_file = filename
    
    grd.verbose = False
    
    grd.read_grid()
    
    for i in [0,1,2,-3,-2,-1]:
        g = grd.g[i]
        print( f"block number: {g.num},lx0={g.lx0},lx1={g.lx1},ly0={g.ly0},ly1={g.ly1},lz0={g.lz0},lz1={g.lz1}" )
    
    blocklist = grd.select_blockgrids([-999,999,-10,0,-999,999], mode='within')
    group_list = grd.group_by_range('xz', block_list=blocklist)

    #print( group_list )
    
    xyz_prb = [76,0,10.4]
    
    print(grd.find_probe_xyz( xyz_prb ))

# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()