# -*- coding: utf-8 -*-
'''
@File    :   block.py
@Time    :   2023/02/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Class block
'''

import warnings
import numpy             as     np
import pandas            as     pd

from   .io_binary        import read_int_bin
from   .io_binary        import read_flt_bin
from   .io_binary        import read_3Dflt_bin

from   .grid             import GridData
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
        can initialized as void or from a binary file
        
        file       : opened file object
        block_list : list of blocks's numbers
        n_var      : number vars in the stored blocks
        vars       : list of selected variable name strings
        var_indx   : list of indexes of selected variables in data chunk
        """
        
        if file is None:
            self.g = GridBlock()
        
        else:
            self._init_from_file( file, block_list, n_var, vars, var_indx, verbose )
            
            # init grid points
            self.g = GridBlock()
    
    def _init_from_file( self, file, block_list, n_var, vars, var_indx, verbose ):
        
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
        self.df_wall = None
        
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


# ----------------------------------------------------------------------
# >>> Drop ghost cells                                             ( 2 )
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

    def drop_ghost( self, buff=3, mode='symmetry' ):
        
        """
        buff: number of ghost layers to be dropped \n
        mode: 'symmetry' or 'oneside'
        """
        
# ----- check if the block is filled

        if self.df.empty:
            warnings.warn(f"drop_ghost: BlockData {self.num} is empty." +
                          f"x,y,z = {self.g.lx0},{self.g.ly0},{self.g.lz0}")
        
# ----- init a cleaned block

        bl_clean = BlockData()
        
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
        
        # decide buffer based on mode
        
        if mode == 'symmetry':  buffl = buff; buffr = buff
        
        elif mode == 'oneside': buffl = buff; buffr = buff-1
        
        # drop ghost cell coordinates
        
        if self.g is not None:
            
            if block_type == 'block':
                bl_clean.g.gx = self.g.gx[buffl:-buffr]
                bl_clean.g.gy = self.g.gy[buffl:-buffr]
                bl_clean.g.gz = self.g.gz[buffl:-buffr]
            elif block_type == 'X':
                bl_clean.g.gy = self.g.gy[buffl:-buffr]
                bl_clean.g.gz = self.g.gz[buffl:-buffr]
            elif block_type == 'Y':
                bl_clean.g.gx = self.g.gx[buffl:-buffr]
                bl_clean.g.gz = self.g.gz[buffl:-buffr]
            elif block_type == 'Z':
                bl_clean.g.gx = self.g.gx[buffl:-buffr]
                bl_clean.g.gy = self.g.gy[buffl:-buffr]

        # drop ghost cells
        # Notice the binary data storage order [n_var, Z, Y, X]
        # Reshape to chunk -> slice -> reshape to vector
        
        if block_type == 'block':
            
            Nx = N1 - (buffl+buffr)
            Ny = N2 - (buffl+buffr)
            Nz = N3 - (buffl+buffr)
            
            sol_buff = (self.df.values).T.reshape( n_var, N3, N2, N1 )
            sol_buff = sol_buff[ :, buffl:-buffr, buffl:-buffr, buffl:-buffr ]

        elif block_type == 'X':
            
            Nx = 1
            Ny = N2 - (buffl+buffr)
            Nz = N3 - (buffl+buffr)
            
            sol_buff = (self.df.values).T.reshape( n_var, N3, N2 )
            sol_buff = sol_buff[ :, buffl:-buffr, buffl:-buffr ]

        elif block_type == 'Y':
            
            Nx = N1 - (buffl+buffr)
            Ny = 1
            Nz = N3 - (buffl+buffr)

            sol_buff = (self.df.values).T.reshape( n_var, N3, N1 )
            sol_buff = sol_buff[ :, buffl:-buffr, buffl:-buffr ]
            
        elif block_type == 'Z':
            Nx = N1 - (buffl+buffr)
            Ny = N2 - (buffl+buffr)
            Nz = 1
        
            sol_buff = (self.df.values).T.reshape( n_var, N2, N1 )
            sol_buff = sol_buff[ :, buffl:-buffr, buffl:-buffr ]
        
        # Reshape the matrix ( no matter 2D or 3D) to long vectors
                
        sol_buff = sol_buff.reshape(( n_var, Nx*Ny*Nz ))
        
        df = pd.DataFrame(sol_buff.T, columns=self.df.columns)
        
        # Append to snap_data_clean
        
        bl_clean.num = self.num
        bl_clean.npx = Nx
        bl_clean.npy = Ny
        bl_clean.npz = Nz
        bl_clean.np  = Nx*Ny*Nz
        bl_clean.df  = df
        
        return bl_clean

# ----------------------------------------------------------------------
# >>> compute gradients                                            ( 3 )
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

    def compute_gradients_block( self, grads:list ):
        
        """
        grads: list of strings, choose from 
        ['grad_rho', 'laplacian', 'grad_p', 'vorticity','Q_cr','lambda2','div',
        'grad_rho_mod','Ducros']
        """
        
        df = self.df
        g  = self.g

        T = np.array( df['T'] )
        p = np.array( df['p'] )
        
        R = 287.0508571
        df['rho'] = p / (T*R)
        
# ----- check if the block is 3D or 2D, if 2D, which type

        if any(x==1 for x in (self.npx,self.npy,self.npz)):
            if self.npx == 1: block_type = 'X'
            if self.npy == 1: block_type = 'Y'
            if self.npz == 1: block_type = 'Z'
        else:
            block_type = 'block'

# ----- compute magnitude of density gradient

        if 'grad_rho' in grads or 'grad_rho_mod' in grads:
            
            rho = np.array( df['rho'] )
            
            if block_type == 'block':
            
                rho = rho.reshape( self.npz, self.npy, self.npx )
                drho_dz  = np.gradient(rho, g.gz, axis=0)
                drho_dy  = np.gradient(rho, g.gy, axis=1)
                drho_dx  = np.gradient(rho, g.gx, axis=2)
                grad_rho = np.sqrt( drho_dx**2 + drho_dy**2 + drho_dz**2 )
            
            elif block_type == 'X':
                
                rho      = rho.reshape( self.npz, self.npy )
                drho_dz  = np.gradient(rho, g.gz, axis=0)
                drho_dy  = np.gradient(rho, g.gy, axis=1)
                grad_rho = np.sqrt( drho_dy**2 + drho_dz**2 )
                
            elif block_type == 'Y':
                
                rho      = rho.reshape( self.npz, self.npx )
                drho_dz  = np.gradient(rho, g.gz, axis=0)
                drho_dx  = np.gradient(rho, g.gx, axis=1)
                grad_rho = np.sqrt( drho_dx**2 + drho_dz**2 )
                
            elif block_type == 'Z':
                
                rho      = rho.reshape( self.npy, self.npx )
                drho_dy  = np.gradient(rho, g.gy, axis=0)
                drho_dx  = np.gradient(rho, g.gx, axis=1)
                grad_rho = np.sqrt( drho_dx**2 + drho_dy**2 )
                
            df['grad_rho'] = grad_rho.flatten()
            
# ----- compute Laplacian of the density

        if 'laplacian' in grads:
            
            if 'grad_rho' not in df.columns:
                raise ValueError("Include 'grad_rho' first!")
            
            if block_type == 'block':
                
                laplacian = np.gradient(drho_dz, g.gz, axis=0) \
                          + np.gradient(drho_dy, g.gy, axis=1) \
                          + np.gradient(drho_dx, g.gx, axis=2)
                df['laplacian'] = laplacian.flatten()
            
            elif block_type == 'X':
                
                laplacian = np.gradient(drho_dz, g.gz, axis=0) \
                          + np.gradient(drho_dy, g.gy, axis=1)
                df['laplacian'] = laplacian.flatten()

            elif block_type == 'Y':
                
                laplacian = np.gradient(drho_dz, g.gz, axis=0) \
                          + np.gradient(drho_dx, g.gx, axis=1)
                df['laplacian'] = laplacian.flatten()
                
            elif block_type == 'Z':
                
                laplacian = np.gradient(drho_dy, g.gy, axis=0) \
                          + np.gradient(drho_dx, g.gx, axis=1)
                df['laplacian'] = laplacian.flatten()


# ----- compute magnitude of pressure gradient

        if 'grad_p' in grads:
            
            p = np.array( df['p'] )
            
            if block_type == 'block':
                
                p      = p.reshape( self.npz, self.npy, self.npx )
                dp_dz  = np.gradient(p, g.gz, axis=0)
                dp_dy  = np.gradient(p, g.gy, axis=1)
                dp_dx  = np.gradient(p, g.gx, axis=2)
                grad_p = np.sqrt( dp_dx**2 + dp_dy**2 + dp_dz**2 )
            
            elif block_type == 'X':
                
                p      = p.reshape( self.npz, self.npy )
                dp_dz  = np.gradient(p, g.gz, axis=0)
                dp_dy  = np.gradient(p, g.gy, axis=1)
                grad_p = np.sqrt( dp_dy**2 + dp_dz**2 )
                
            elif block_type == 'Y':
                
                p      = p.reshape( self.npz, self.npx )
                dp_dz  = np.gradient(p, g.gz, axis=0)
                dp_dx  = np.gradient(p, g.gx, axis=1)
                grad_p = np.sqrt( dp_dx**2 + dp_dz**2 )
                
            elif block_type == 'Z':
                
                p      = p.reshape( self.npy, self.npx )
                dp_dy  = np.gradient(p, g.gy, axis=0)
                dp_dx  = np.gradient(p, g.gx, axis=1)
                grad_p = np.sqrt( dp_dx**2 + dp_dy**2 )
                
            df['grad_p'] = grad_p.flatten()
            
# -- if vorticity, div, or Ducros are selected, compute velocity gradient first

        gradV = ['vorticity','Q_cr','lambda2','div','grad_rho_mod','Ducros']
        if any(_ in grads for _ in gradV):

            if block_type == 'block':
                
                u = np.array( df['u'] ).reshape( self.npz, self.npy, self.npx )
                v = np.array( df['v'] ).reshape( self.npz, self.npy, self.npx )
                w = np.array( df['w'] ).reshape( self.npz, self.npy, self.npx )
                
                du_dx = np.gradient( u, g.gx, axis=2 )
                du_dy = np.gradient( u, g.gy, axis=1 )
                du_dz = np.gradient( u, g.gz, axis=0 )
                dv_dx = np.gradient( v, g.gx, axis=2 )
                dv_dy = np.gradient( v, g.gy, axis=1 )
                dv_dz = np.gradient( v, g.gz, axis=0 )
                dw_dx = np.gradient( w, g.gx, axis=2 )
                dw_dy = np.gradient( w, g.gy, axis=1 )
                dw_dz = np.gradient( w, g.gz, axis=0 )
            
            elif block_type == 'X':
                
                v = np.array( df['v'] ).reshape( self.npz, self.npy )
                w = np.array( df['w'] ).reshape( self.npz, self.npy )
                
                dv_dy = np.gradient( v, g.gy, axis=1 )
                dv_dz = np.gradient( v, g.gz, axis=0 )
                dw_dy = np.gradient( w, g.gy, axis=1 )
                dw_dz = np.gradient( w, g.gz, axis=0 )

            elif block_type == 'Y':
                
                u = np.array( df['u'] ).reshape( self.npz, self.npx )
                w = np.array( df['w'] ).reshape( self.npz, self.npx )
                
                du_dx = np.gradient( u, g.gx, axis=1 )
                du_dz = np.gradient( u, g.gz, axis=0 )
                dw_dx = np.gradient( w, g.gx, axis=1 )
                dw_dz = np.gradient( w, g.gz, axis=0 )
            
            elif block_type == 'Z':
                
                u = np.array( df['u'] ).reshape( self.npy, self.npx )
                v = np.array( df['v'] ).reshape( self.npy, self.npx )
                
                du_dx = np.gradient( u, g.gx, axis=1 )
                du_dy = np.gradient( u, g.gy, axis=0 )
                dv_dx = np.gradient( v, g.gx, axis=1 )
                dv_dy = np.gradient( v, g.gy, axis=0 )
            
# ----- compute vorticity

        if 'vorticity' in grads:
            
            if block_type == 'block':
                
                w1 = dw_dy - dv_dz
                w2 = du_dz - dw_dx
                w3 = dv_dx - du_dy
                
                df['w1'] = w1.flatten()
                df['w2'] = w2.flatten()
                df['w3'] = w3.flatten()
                df['vorticity'] = np.sqrt( w1**2 + w2**2 + w3**2 ).flatten()
            
            if block_type == 'X':
                w1 = dw_dy - dv_dz
                df['w1'] = w1.flatten()
            
            if block_type == 'Y':
                w2 = du_dz - dw_dx
                df['w2'] = w2.flatten()
            
            if block_type == 'Z':
                w3 = dv_dx - du_dy
                df['w3'] = w3.flatten()
                
# ----- compute Q-cr

        if 'Q_cr' in grads:
            
            if block_type != 'block':
                raise ValueError("Q_cr is only for block type!")
            
            t_comp = 1/3*(du_dx + dv_dy + dw_dz)

            Q_cr   = -0.5*( (du_dx-t_comp)**2 +
                            (dv_dy-t_comp)**2 +
                            (dw_dz-t_comp)**2 ) \
                     - (du_dy*dv_dx + du_dz*dw_dx + dv_dz*dw_dy)
            
            df['Q_cr'] = Q_cr.flatten()  

# ----- compute lambda2

        if 'lambda2' in grads:
            
            if block_type != 'block':
                raise ValueError("lambda2 is only for block type!")
            
            J = np.array([[du_dx,du_dy,du_dz],
                          [dv_dx,dv_dy,dv_dz],
                          [dw_dx,dw_dy,dw_dz]])
            
            S     = 0.5 * ( J + np.transpose(J,(1,0,2,3,4)) )
            Omega = 0.5 * ( J - np.transpose(J,(1,0,2,3,4)) )
            
            eigs = np.zeros( (3,self.npz,self.npy,self.npx) )

            for i in range(self.npz):
                for j in range(self.npy):
                    for k in range(self.npx):
                        eigs_unsorted = np.linalg.eigvals(S[:,:,i,j,k]**2+Omega[:,:,i,j,k]**2)
                        eigs[:,i,j,k] = np.sort(eigs_unsorted)[::-1]
                        
            df['lambda2'] = eigs[1,:,:,:].flatten()   

# ----- compute the divergence of velocity

        if 'div' in grads:
            
            if block_type != 'block':
                raise ValueError("div is only for block type!")

            div = du_dx + dv_dy + dw_dz
            df['div'] = div.flatten()

# ----- compute modified grad_rho

        if 'grad_rho_mod' in grads:

            if block_type != 'block':
                raise ValueError("grad_rho_mod is only for block type!")
            
            grad_rho = np.array( df['grad_rho'] )
            div      = np.array( df['div']      )
                        
            # shock region, div should be negative, and vorticity should be quite small;
            # in boundary layer region, div fluctuates but vorticity is relatively large.
            # so, we can use 'sensor' to filter out the boundary layer region.
            
            sensor = div / ( np.array(df['vorticity']) + 1e-10)

            grad_rho[np.where( sensor > -0.8 )] = 0.0

            df['sensor']       = sensor
            df['grad_rho_mod'] = grad_rho
            
# ----- compute Ducros sensor

        if 'Ducros' in grads:

            if block_type != 'block':
                raise ValueError("Ducros is only for block type!")
            
            div     = du_dx + dv_dy + dw_dz
            curl_sq = (dw_dy-dv_dz)**2 + (du_dz-dw_dx)**2 + (dv_dx - du_dy)**2
            Ducros  = div**2 / ( div**2 + curl_sq**2 + 1e-30)

            df['Ducros'] = Ducros.flatten()

# ----- update dataframe

        self.df    = df
        self.n_var = len( df.columns )
        
        
# ----------------------------------------------------------------------
# >>> Define a subclass SnapBlock                                  ( 4 )
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
            self._init_from_file( file, block_list, n_var, vars, snap_with_gx,
                                 type, kind)
        

    def _init_from_file( self, file, block_list, var_read, n_var, vars, snap_with_gx, 
                        type, kind=4):
        
        """
        After initialization, self.g only contains cell center coordinates.
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
        self.n_var = n_var
        
        # the variables that will be read; if var_read is None, read all
        
        if var_read is None: var_read = vars
        
        # empty grids points list
        
        n_grid = 0
        self.g = GridBlock()
        
        # empty dataframe for datachunk

        self.df = pd.DataFrame(columns=var_read)
        
        # matrix of friction projection on x-z plane

        self.df_fric = None
        self.df_wall = None
        
# ----- read global block number and block dimensions
        
        self.num = read_int_bin( file.read(self.sin), self.sin )
        
        if self.num == 0:
            raise ValueError(f"Found a block with num==0 at {pos_start}.")
        
        dim      = read_int_bin( file.read(3*self.sin), self.sin )
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
            pos      += self.npx*kind
            
            self.g.gy = read_3Dflt_bin( pos, file, 1, self.npy, 1, kind )
            pos      += self.npy*kind
            
            self.g.gz = read_3Dflt_bin( pos, file, 1, 1, self.npz, kind )
            pos      += self.npz*kind
            
            n_grid    = self.npx + self.npy + self.npz
        
        if snap_with_gx and type == 'slice':
            
            # slic_type == 'X'
            
            if self.npx == 1: 
                
                self.g.gx = None
                
                self.g.gy = read_3Dflt_bin( pos, file, 1,self.npy,1, kind )
                pos      += self.npy*kind

                self.g.gz = read_3Dflt_bin( pos, file, 1,1,self.npz, kind )
                pos      += self.npz*kind
                
                n_grid    = self.npy + self.npz

            # slic_type == 'Y' or 'W'
            
            elif self.npy == 1:

                self.g.gx = read_3Dflt_bin( pos, file, self.npx,1,1, kind )
                pos      += self.npx*kind
                
                self.g.gy = None

                self.g.gz = read_3Dflt_bin( pos, file, 1,1,self.npz, kind )
                pos      += self.npz*kind
                
                n_grid    = self.npx + self.npz
            
            # slic_type == 'Z'
            
            elif self.npz == 1:
                
                self.g.gx = read_3Dflt_bin( pos, file, self.npx,1,1,kind )
                pos      += self.npx*kind

                self.g.gy = read_3Dflt_bin( pos, file, 1,self.npy,1,kind )
                pos      += self.npy*kind
                
                self.g.gz = None
                
                n_grid    = self.npx + self.npy

# ----- record the starting position of variable chunk
            
        self.pos_var_start = pos
                
# ----- read variable data chunk

        # if this block chunk will be read?
        # by default, to_fill is False
        self.to_fill = False
        
        if (block_list is None) or (self.num in block_list):
            self.to_fill = True

        # primitive variables ['u', 'v', 'w', 'T', 'p'] with current setting
        # only read variables in var_read
        if self.to_fill:
            
            for var in vars:
            
                if var in var_read:
                    buff         = read_flt_bin( file.read(self.np*kind), kind )
                    self.df[var] = buff
                
                else:
                    file.seek( self.np*kind, 1 )
                    
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
        self.df_wall = None
        
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


# ----------------------------------------------------------------------
# >>> Class SolutionBlock                                         (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/05/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

class SolutionBlock(BlockData):
    
    # size of data type
    sin, sfl, slg = 4, 8, 4
    
    def __init__( self, fs=None, block_list=None, n_vars=None, vars=None, 
                        G=None, verbose=False ):
        
        if fs is None:
            self.g = GridBlock()
        
        else: 
            self._init_from_file( fs, G, block_list, n_vars, vars, verbose)
    
    def _init_from_file( self, fs, G:GridData, block_list:list, n_vars, 
                               vars:list, verbose):

        # start position of a file pointer
        pos_start = fs.tell()
        
        # empty list for future use
        self.df = pd.DataFrame(columns=vars)
        
        # matrix of friction projection on x-z plane
        self.df_fric = None
        self.df_wall = None
        
        # read global block number(index)
        
        self.num = read_int_bin( fs.read(self.sin), self.sin )
        print(self.num)
        
        self.g   = G.g[self.num-1]
        
        self.npx = self.g.nx + 6
        self.npy = self.g.ny + 6
        self.npz = self.g.nz + 6
        self.np  = self.g.np
        
        # if this block will be read? by default, to_fill is False
        self.to_fill = False
        if self.num in block_list: self.to_fill = True
        
        if verbose:
            print(f"read in block {self.num},",end='')
            print(f"dimension {self.npx} {self.npy} {self.npz}. {self.to_fill}")
        
        # read in data into dataframe
        if self.to_fill:
            
            for var in vars:
                
                pos  = pos_start + self.sin + vars.index(var)*self.np*self.sfl
                fs.seek( pos )
                
                buff = read_flt_bin( fs.read(self.np*self.sfl), self.sfl )
                self.df[var] = buff
        
        # calculate the block data size in byte
        
        self.size = self.np*n_vars*self.sfl + self.sin
        
        # move file pointer to the end of the current block
        
        fs.seek( pos_start + self.size )

