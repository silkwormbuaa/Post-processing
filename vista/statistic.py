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
import gc
import sys
import pickle
import numpy             as     np
import pandas            as     pd
from   copy              import deepcopy

from   .grid             import GridData
from   .block            import BlockData
from   .tools            import find_indices
from   .io_vtk           import write_vtm_file
from   .io_vtk           import create_multiblock_dataset
from   .io_vtk           import add_var_vtkRectilinearGrid
from   .io_vtk           import create_3d_vtkRectilinearGrid
from   .io_binary        import read_int_bin
from   .io_binary        import read_flt_bin
from   .io_binary        import read_log_bin
from   .io_binary        import write_flt_bin
from   .io_binary        import write_int_bin
from   .io_binary        import write_log_bin

from   .lib.form         import phy
from   .lib.form         import mth

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

    def __init__( self, stat_file=None ):
        
        """
        stat_file: path to statistics.bin.
                   (optional) if None, generate a void Statistic.\n
        
        return     :  self.fsize \n  
        initialize :  self.pos, self.n_var, self.bl, self.verbose
        """
        
        if stat_file is not None:
        
            # directory to the statistic file
            self.file = stat_file

            # file size
            self.fsize = os.stat(self.file).st_size
        
        else: 
            
            self.file  = None
            self.fsize = 0
        
        # file pointer position
        self.pos = 0

        # number of variables
        self.n_var = 0
        
        # type of statistic, 'block' or 'X', 'Y', 'Z'
        self.stat_type = None
        
        # list of blocks
        self.bl       = []
        self.bl_clean = []

        # list of block numbers
        self.bl_nums       = []
        self.bl_nums_clean = []
        
        # Verbose ? 
        self.verbose = False
        
        # Grid3d (including all blocks grids from inca_grid.bin file)
        self.grid3d        = None
        
        # Dataframe for friction and pressure projection 
        self.df_fric = None
        self.df_wall = None
        

# ----------------------------------------------------------------------
# >>> read statistic                                          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_statistic( self, block_list, vars_in ):
        
        """
        block_list : list of selected blocks' numbers \n
        vars_in    : list of variable names \n
        """
        
        with open( self.file, 'rb' ) as f:
            
            self.read_stat_header( f )
            self.read_stat_body( f, block_list, vars_in )

        # determined stat_type
        
        bl = self.bl[ 0 ]
        
        if   bl.npx == 1: self.stat_type = 'X'
        elif bl.npy == 1: self.stat_type = 'Y'
        elif bl.npz == 1: self.stat_type = 'Z'
        else            : self.stat_type = 'block'
        
        print(f"StatisticData {self.file} (type: {self.stat_type}) is read.\n")
        sys.stdout.flush()
        
        
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
        
        """
        file:  opened file object\n
        
        return: \n
            file header info ( self.verbose=True to show) and
            self.pos will be updated
        """

        # Internal variables

        sfl = 8
        sin = 4
        slg = 4

        # Read file format

        self.format = read_int_bin(file.read(sin),sin)

        self.pos += sin


        if self.format == 1:

            # Read number of samples, sample step & start step

            self.n_samples = read_int_bin(file.read(sin),sin)
            self.sample_step = read_int_bin(file.read(sin),sin)
            self.start_step  = read_int_bin(file.read(sin),sin)

            self.pos += int(3*sin)

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
            
        
        elif self.format == 2:

            # Read number of samples, sample step & start step

            self.n_samples   = read_int_bin(file.read(sin),sin)
            self.sample_step = read_int_bin(file.read(sin),sin)
            self.start_step  = read_int_bin(file.read(sin),sin)

            self.pos += int(3*sin)
            
            # Read number of transported variables
            
            self.npv      = read_int_bin(file.read(sin),sin) 
            self.nscalars = read_int_bin(file.read(sin),sin)
            self.nspecies = read_int_bin(file.read(sin),sin)
            
            self.pos += int(3*sin)

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

        else:
            raise ValueError('Header format in statistics.bin not supported')

        # if want to see header, set verbose True.            
        if (self.verbose is True) : 
            print('format           '+ str(self.format)           + '\n')
            print('n_samples        '+ str(self.n_samples)        + '\n')
            print('start_step       '+ str(self.start_step)       + '\n')
            print('sample_step      '+ str(self.sample_step)      + '\n')
            print('sample_time      '+ str(self.sample_time)      + '\n')
            print('start_time       '+ str(self.start_time)       + '\n')
            print('meanvalues       '+ str(self.meanvalues)       + '\n')
            print('doublecorr       '+ str(self.doublecorr)       + '\n')
            print('triplecorr       '+ str(self.triplecorr)       + '\n')
            print('quadruplecorr    '+ str(self.quadruplecorr)    + '\n')
            print('autocorr         '+ str(self.autocorr)         + '\n')
            print('mean_invar       '+ str(self.mean_invar)       + '\n')
            print('schlieren        '+ str(self.schlieren)        + '\n')
            print('cavitation_stats '+ str(self.cavitation_stats) + '\n')
            print('vapor_gas_stats  '+ str(self.vapor_gas_stats)  + '\n')
            print('rste             '+ str(self.rste)             + '\n')
            print('thermo           '+ str(self.thermo)           + '\n')
            print('visc_diff        '+ str(self.visc_diff)        + '\n')
            print('header_size      '+ str(self.header_size)      + '\n')


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

    def read_stat_body( self, file , bl_list, sel_vars ):
        
        """
        file       :  opened file object \n
        bl_list    :  list of selected blocks' numbers \n 
        vars       :  list of variable name strings \n
        
        return     :  self.bl (list of BlockData)
        """
        verbose = self.verbose
        end_of_file = False
               
        # New position after reading headers
        
        file.seek( self.header_size )
        
        self.pos = self.header_size
        
        # vars determine the index of data sequences to read
        
        vars_indx = self.vars_to_indx( sel_vars )
        
        # total number of variables in the block
        
        n_var = self.n_var
        
        while not end_of_file:
        
            # read in block one by one
            # only blocks in fill list will be filled with data chunk
        
            self.bl.append(BlockData(file, bl_list, n_var, sel_vars, vars_indx, verbose))
            
            # collect the bl_nums and match grid to block
            
            self.bl_nums.append(self.bl[-1].num)
            
            if self.grid3d is not None:
                self.bl[-1].g = self.grid3d.g[self.bl[-1].num-1]
            
            self.pos = self.pos + self.bl[-1].size

            if self.pos >= self.fsize: end_of_file = True
            
        # get list of block numbers
        
        self.bl_nums = [block.num for block in self.bl]


# ----------------------------------------------------------------------
# >>> vars to indx                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/01  - created
#
# Desc
#   - from variables list to sequence index list within in BlockData
# ----------------------------------------------------------------------
    
    def vars_to_indx( self, vars ):
        
        """
        vars: list of variable name strings
        """
        
        # list of variables for inquiry
        
        mean_ls = ['u','v','w','rho','rhoE','p','T','mu']
        
        cor2_ls = ['uu','uv','uw','urho'  ,'urhoE'   ,'up'   ,'uT'   ,'umu'   ,
                        'vv','vw','vrho'  ,'vrhoE'   ,'vp'   ,'vT'   ,'vmu'   ,
                             'ww','wrho'  ,'wrhoE'   ,'wp'   ,'wT'   ,'wmu'   ,
                                  'rhorho','rhorhoE' ,'rhop' ,'rhoT' ,'rhomu' ,
                                           'rhoErhoE','rhoEp','rhoET','rhoEmu',
                                                      'pp'   ,'pT'   ,'pmu'   ,
                                                              'TT'   ,'Tmu'   ,
                                                                      'mumu'  ]   
        
        # index list of vars
        
        indx = []
        
        # from variables list to corresponding indexes list
        
        for var in vars:
            
            displace = 0
            
            # mean variables
            
            if var in mean_ls:
                indx.append( displace + mean_ls.index(var) )

            if self.meanvalues: displace += 8
            
            # double correlated variables
            
            if var in cor2_ls:
                indx.append( displace + cor2_ls.index(var) )
                
            if self.doublecorr: displace += 36
        
        return indx        


# ----------------------------------------------------------------------
# >>> Match grid                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/01  - created
#
# Desc
#   - match data with grids
# ----------------------------------------------------------------------
    
    def match_grid( self, block_list, G, add_to_df=True ):
        
        """
        bl.g => g. Necessary if grid3d was not set before read_stat()
        
        block_list : list of selected blocks' numbers \n
        G          : GridData object \n
        add_to_df  : add x,y,z and vol_fra to self.bl[].df, if not, just connect
                     g to self.bl[index]\n
        
        return : x,y,z and vol_fra are added to self.bl[].df \n
        """
        
        stat_type = self.stat_type
        
        for num in block_list:
            
            g = G.g[num-1]
            
            bl = self.bl[self.bl_nums.index(num)]
            
            bl.g = g
            
            if add_to_df:
                if stat_type == 'block':
                
                    X, Y, Z  = np.meshgrid( g.gx, g.gy, g.gz, indexing='ij' )
                    hx,hy,hz = np.meshgrid( g.hx, g.hy, g.hz, indexing='ij' )
                    
                    bl.df['x']  = X.T.flatten()
                    bl.df['y']  = Y.T.flatten()
                    bl.df['z']  = Z.T.flatten()
                    bl.df['hx'] = hx.T.flatten()
                    bl.df['hy'] = hy.T.flatten()
                    bl.df['hz'] = hz.T.flatten()
                    
                    
                    # adding vol_fra !! original vol_fra has i,j,k order, 
                    # should be transpose as k,j,i
                    if g.vol_fra is not None:
                        bl.df['vol_fra'] = np.ravel( g.vol_fra.T )
                
                    if g.vol is not None:
                        bl.df['vol'] = np.ravel( g.vol )

                elif stat_type == 'X':
                    
                    Z, Y  = np.meshgrid( g.gz, g.gy, indexing='ij' )
                    hz,hy = np.meshgrid( g.hz, g.hy, indexing='ij' )

                    bl.df['y']  = Y.flatten()
                    bl.df['z']  = Z.flatten()
                    bl.df['hy'] = hy.flatten()
                    bl.df['hz'] = hz.flatten()
                
                elif stat_type == 'Y':
                    
                    Z, X  = np.meshgrid( g.gz, g.gx, indexing='ij' )
                    hz,hx = np.meshgrid( g.hz, g.hx, indexing='ij' )
                    
                    bl.df['x']  = X.flatten()
                    bl.df['z']  = Z.flatten()
                    bl.df['hx'] = hx.flatten()
                    bl.df['hz'] = hz.flatten()
                
                elif stat_type == 'Z':
                    
                    Y, X  = np.meshgrid( g.gy, g.gx, indexing='ij' )
                    hy,hx = np.meshgrid( g.hy, g.hx, indexing='ij' )
                    
                    bl.df['x']  = X.flatten()
                    bl.df['y']  = Y.flatten()
                    bl.df['hx'] = hx.flatten()
                    bl.df['hy'] = hy.flatten()
                    

                    
# ----------------------------------------------------------------------
# >>> drop ghost cells                                          (Nr.)
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

    def drop_ghost( self, block_list, buff=3, mode='symmetry' ):
        
        """
        block_list : list of blocks that are going to drop ghost cells \n
        buff       : number of buffer layers \n
        mode       : 'symmetry' or 'oneside' \n
        return     : self.bl_clean[]
        """
# ----- check if data is available

        if len(self.bl) == 0:
            raise ValueError('No data in the statistics.')

# ----- drop ghost cells and store data in the self.bl_clean[]

        else:
            
            del self.bl_clean
            gc.collect()
            
            self.bl_clean = []
            
            for block in self.bl:
                
                if block.num not in block_list:
                    continue
                
                self.bl_clean.append( block.drop_ghost(buff=buff,mode=mode) )

            self.bl_nums_clean = [block.num for block in self.bl_clean]

            print('Ghost cells in selected blocks are dropped.')
        

# ----------------------------------------------------------------------
# >>> compute Mach number                                         (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/10  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_vars( self, block_list:list, vars_new:list ):
        
        """
        block_list : list of blocks that are going to compute new variables\n
        vars_new   : list of str representing new vars 
                     ['mach','RS','p`','mu','favre_velocity']\n
        
        return : corresponding variables are added to self.bl[].df
        """
        
        for num in block_list:
            
            df = self.bl[self.bl_nums.index(num)].df
            
# --------- compute Mach number 

            if "mach" in vars_new:
                
                u = np.array( df['u'] )
                v = np.array( df['v'] )
                w = np.array( df['w'] )
                T = np.array( df['T'] )
                gamma = 1.40
                R = 287.0508571
                
                mach = np.sqrt( u*u+v*v+w*w ) / np.sqrt( gamma*R*T )
                
                self.bl[self.bl_nums.index(num)].df['mach'] = mach

# ---------- compute turbulent kinetic energy

            if "RS" in vars_new:
                
                uu  = np.array(df['uu']) - np.array(df['u'])*np.array(df['u'])
                vv  = np.array(df['vv']) - np.array(df['v'])*np.array(df['v'])
                ww  = np.array(df['ww']) - np.array(df['w'])*np.array(df['w'])
                uv  = np.array(df['uv']) - np.array(df['u'])*np.array(df['v'])
                tke = uu + vv + ww
                
                df['u`u`'] = uu 
                df['v`v`'] = vv
                df['w`w`'] = ww
                df['u`v`'] = uv
                df['tke']  = tke  

# ---------- compute pressure fluctuation

            if "p`" in vars_new:
                
                # in solid cells, may appear negative values, so use abs.
                p_fluc = abs(np.array(df['pp']) - np.array(df['p'])**2)
                df['p`'] = np.sqrt( p_fluc )

# ---------- compute viscosity mu

            if "mu" in vars_new:
                
                mu = np.array( df['T'].apply(phy.sutherland) )
                df['mu'] = mu

# ---------- compute favre average velocity

            if "favre_velocity" in vars_new:
                
                rho           = np.array( df['rho']  )
                df['rho']     = np.array( df['rho']  )
                df['u_favre'] = np.array( df['urho'] ) / rho
                df['v_favre'] = np.array( df['vrho'] ) / rho
                df['w_favre'] = np.array( df['wrho'] ) / rho
                

# ----------------------------------------------------------------------
# >>> compute gradients of variables                              (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/11  - created
# 2024/02/22  - computation method is within Block
#
# Desc
#   - need to be improved when dealing with near-wall cells
# ----------------------------------------------------------------------
                
    def compute_gradients( self, block_list:list, grads:list ):
        
        """
        block_list: list of blocks that are going to compute gradients\n
        grads: list of str representing gradients from
        ['grad_rho', 'laplacian', 'grad_p', 'vorticity','Q_cr','lambda2','div']\n
        """
        
        for num in block_list:
            
            block = self.bl[self.bl_nums.index(num)]
            
            block.compute_gradients_block( grads )
        
        print(f"Gradients {grads} are computed.\n")


# ----------------------------------------------------------------------
# >>> compute source term of secondary flow                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_source_terms( self, block_list:list, G:GridData, buff=3):
        
        """
        block_list: list of blocks that are going to compute gradients\n
        G: GridData instance of corresponding case\n
        
        return: self.bl[num-1].df['S1'] (...['S2'])
        """
        
        for num in block_list:
            
            df = self.bl[self.bl_nums.index(num)].df
            g  = G.g[num-1]

            vars_exist = df.columns
            
            if 'v`v`' not in vars_exist:
                vv = np.array(df['vv']) - np.array(df['v'])*np.array(df['v'])
            else: vv = np.array( df['v`v`'] )
            
            if 'w`w`' not in vars_exist:
                ww = np.array(df['ww']) - np.array(df['w'])*np.array(df['w'])
            else: ww = np.array( df['w`w`'] )
            
            vw = np.array(df['vw']) - np.array(df['v'])*np.array(df['w'])
            
            npx = g.nx + buff*2
            npy = g.ny + buff*2
            npz = g.nz + buff*2
            
# --------- compute S1 = dd(v'v'-w'w')/dydz

            temp = np.array( vv-ww ).reshape( npz, npy, npx )
            
            S1 = np.gradient( np.gradient(temp, g.gy, axis=1), g.gz, axis=0)

            df['S1'] = S1.flatten()
            
# --------- compute S2 = dd(v'w')/dzdz - dd(v'w')/dydy

            vw = vw.reshape( npz, npy, npx)
            
            S2 = np.gradient( np.gradient(vw, g.gz, axis=0), g.gz, axis=0) \
               - np.gradient( np.gradient(vw, g.gy, axis=1), g.gy, axis=1)
            
            df['S2'] = S2.flatten()
            
            df['S']  = S1.flatten() + S2.flatten()
            

# ----------------------------------------------------------------------
# >>> compute profile                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_profile( self, block_list, bbox, vars:list, 
                         RS=True, outfile=False, roughwall=True ):
        
        """
        block list: selected block list (within bounding box)\n
        bbox: bounding box coordinates [xmin,xmax,ymin,ymax,zmin,zmax]\n
        vars: variable names list\n
        RS: set True to compute Reynolds Stresses, otherwise drop them\n
        outfile: assign outfile name or use default 'profile_spanwisemean.dat'
        roughwall: bool value, if True, set all values at wall to zero
        """
        
# ----- collect data frame from all filled blocks
        
        df = pd.concat( [self.bl_clean[self.bl_nums_clean.index(num)].df for num in block_list] )

        # reset indexes in case repeated indexes from different blocks
        
        df.reset_index( drop=True, inplace=True )
        
        print(df)
        
        
# ----- slim down data with bbox
        
        df.drop( df[ (df['x'] < bbox[0]) | (df['x'] > bbox[1]) |
                     (df['y'] < bbox[2]) | (df['y'] > bbox[3]) |
                     (df['z'] < bbox[4]) | (df['z'] > bbox[5]) ].index,
                 inplace=True )
        
        print(df)
        
# ----- compute turbulence (Reynolds Stress) and turbulent kinetic energy
        
        if RS:
            df['u`u`']= np.array(df['uu']) - np.array(df['u'])**2
            df['v`v`']= np.array(df['vv']) - np.array(df['v'])**2
            df['w`w`']= np.array(df['ww']) - np.array(df['w'])**2
            df['u`v`']= np.array(df['uv']) - np.array(df['u'])*np.array(df['v'])
            df['tke'] = df['u`u`'] + df['v`v`'] + df['w`w`']
            
            # delete some var in vars and add some
            
            vars =[var for var in vars if var not in ['uu','vv','ww','uv','pp']]
            vars += ['u`u`','v`v`','w`w`','u`v`','tke']
            
            print(df)

# ----- compute total temperature and total pressure

        ke = 0.5*(np.array(df['uu']) + np.array(df['vv']) + np.array(df['ww']))
        R  = 287.0508571
        gamma = 1.4
        Cp = R*gamma/(gamma-1)
        df['Tt'] = np.array(df['T']) + ke/Cp
        df['pt'] = np.array(df['p'])*(np.array(df['Tt'])/np.array(df['T']))**(gamma/(gamma-1))
        vars += ['Tt','pt']

# ----- do averaging in x-z plane

        # check if 'vol_fra' is in df.columns (so that smooth wall case
        # can be handled )
        if 'vol_fra' not in df.columns:
            df['vol_fra'] = 1.0
        
        ys = np.sort( np.unique( np.array( df['y'] )))
        
        data_chunk = None
        
        vars = ['hy'] + vars
        for y in ys:
            
            df_temp = df[ df['y']==y ]
            
            # only applicable to uniform grid in x-z plane, otherwise need grid volumes
            
            vol     = np.array( df_temp['vol'] )
            vol_fra = np.array( df_temp['vol_fra'] )
            vol_total = np.sum( vol_fra*vol )
            if vol_total < 0.0000001: vol_total = float('inf')
            
            buff = [y]
            
            for var in vars:
                v = np.sum( np.array(df_temp[var])*vol_fra*vol ) / vol_total
                buff.append(v)
                
                # integral of volume fraction in instrinsic averaging
                
                # if abs(y) < 0.2 and var == 'u':
                #     totalval = np.sum( np.array(df_temp[var])*vol_fra )
                #     print(f"{y:12.4f}, {totalval:12.4f}, {np.sum(vol_fra):12.4f},{totalval/np.sum(vol_fra):12.4f}")
            
            if data_chunk is None: data_chunk = [buff]
            else: data_chunk.append( buff )
            
        
        data_chunk = np.array(data_chunk).reshape( len(ys), len(vars)+1 )
        
        df_profile = pd.DataFrame(data_chunk,columns=['y']+vars)
        
        pd.set_option('display.max_rows', None)  # 显示所有行
        
# ----- drop points below wall and set values(except p,rho,T) at wall manually
        
        df_profile.drop( df_profile[ df_profile['p']==0.0 ].index, inplace=True)
        
        df_profile.reset_index( drop=True, inplace=True )
        
        if roughwall:
            for var in vars:
                if var not in ['p','rho','mu','T','Tt','pt']:
                    df_profile.loc[0,var] = 0.0
                
        print("points below wall are droppped and u,v,w are set to zero at wall.")
        print( df_profile )
        
        pd.reset_option('display.max_rows')  # rest max_rows displayed
        
# ----- output profile into txt

        if outfile is False: outfile = 'profile_spanwisemean.dat'

        df_profile.to_string( outfile,  
                              index=False,
                              float_format='%15.7f',
                              justify='left')
         
            
# ----------------------------------------------------------------------
# >>> get slice dataframe (statistics.bin)                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/03  - created
#
# Desc
#   - slice 3D statistic data
#
# ----------------------------------------------------------------------

    def get_slice_df( self, block_list, G, indx_slic, slic_type, buff=3 ):
        
        """
        block_list : list of sliced blocks \n
        G          : corresponding GridData instance \n
        indx_slic  : list of slice indexes on each block \n
        slic_type  : 'X','Y' or 'Z' \n
        
        return : dataframe of sliced results
        """
        
# ----- should somehow check if Statistic is 3D or 2D data?
        
        # pass
        
# ----- do slicing on each bl.df

        for i, num in enumerate( block_list ):
            
            bl   = self.bl[self.bl_nums.index(num)]
            g    = G.g[num-1]
            indx = indx_slic[i]
            
            npx = g.nx + buff*2
            npy = g.ny + buff*2
            npz = g.nz + buff*2
            
            vars = bl.df.columns
            
            data_chunk = None
            
            for var in vars:
                
                data = np.array( bl.df[var] )
                data = data.reshape( npz, npy, npx )
                
                if   slic_type == 'X': data = data[:,:,indx].flatten()
                elif slic_type == 'Y': data = data[:,indx,:].flatten()
                elif slic_type == 'Z': data = data[indx,:,:].flatten()
                
                if data_chunk is None: data_chunk = [data]
                else: data_chunk.append( data )
            
            data_chunk = np.array(data_chunk).T
            
            bl.df = pd.DataFrame(data_chunk,columns=vars)
    
# ----- match grids with data

            if slic_type == 'X':
                
                Y,Z = np.meshgrid( g.gy, g.gz )
                bl.df['y'] = Y.flatten()
                bl.df['z'] = Z.flatten()

            if slic_type == 'Y':
                
                X,Z = np.meshgrid( g.gx, g.gz )
                bl.df['x'] = X.flatten()
                bl.df['z'] = Z.flatten()
                
            if slic_type == 'Z':
                
                X,Y = np.meshgrid( g.gx, g.gy )
                bl.df['x'] = X.flatten()
                bl.df['y'] = Y.flatten()                

# ----- drop ghost cells

            vars = bl.df.columns
            
            data_chunk = None
            
            for var in vars:
                
                data = np.array( bl.df[var] )
                
                if   slic_type == 'X': data = data.reshape( npz, npy )
                elif slic_type == 'Y': data = data.reshape( npz, npx )
                elif slic_type == 'Z': data = data.reshape( npy, npx )
                    
                data = data[buff:-buff,buff:-buff].flatten()

                if data_chunk is None: data_chunk = [data]
                else: data_chunk.append(data)
        
            data_chunk = np.array(data_chunk).T

            bl.df = pd.DataFrame(data_chunk,columns=vars)

# ----- concatenate all dataframes in all selected blocks

        df_slice = pd.concat( [self.bl[self.bl_nums.index(num)].df for num in block_list] )
        
        df_slice.reset_index( drop=True, inplace=True)
        
        return df_slice


# ----------------------------------------------------------------------
# >>> Get probed data                                            (Nr.)
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

        for i, num in enumerate( bl_list ):
            
            bl    = self.bl[self.bl_nums.index(num)]
            g     = G.g[num-1]
            indx  = indx_probed[i]
            
            npx   = g.nx + buff*2
            npy   = g.ny + buff*2
            npz   = g.nz + buff*2
            
            vars = bl.df.columns
            
            data_chunk = None
            
            for var in vars:
                
                data = np.array( bl.df[var] )
                data = data.reshape( npz, npy, npx )
                
                if   probe_type == 'X': data = data[indx[0],indx[1],:]
                elif probe_type == 'Y': data = data[indx[1],:,indx[0]]
                elif probe_type == 'Z': data = data[:,indx[1],indx[0]]
                
                if data_chunk is None: data_chunk = [data]
                else: data_chunk.append(data)
                
            data_chunk = np.array(data_chunk).T
            
            bl.df_probe = pd.DataFrame(data_chunk,columns=vars)
            
# ------ match grids with data

            if probe_type == 'X': bl.df_probe['x'] = g.gx
            if probe_type == 'Y': bl.df_probe['y'] = g.gy
            if probe_type == 'Z': bl.df_probe['z'] = g.gz
        
# ------ drop ghost cells

            bl.df_probe = bl.df_probe.iloc[buff:-buff]
        
# ------ concatenate all dataframe in all selected blocks

        df_probe = pd.concat( [self.bl[self.bl_nums.index(num)].df_probe for num in bl_list] )
        
        df_probe.sort_values(by=[probe_type.lower()],inplace=True)
        
        df_probe.reset_index( drop=True, inplace=True )
        
        return df_probe
        

# ----------------------------------------------------------------------
# >>> wall variables projection                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/23  - created
#
# Desc
#    - compute the friction distribution projected on x-z plane
#    - return a dataframe and output file
# ----------------------------------------------------------------------

    def friction_projection( self,         block_list:list,
                             G:GridData,   cc_df:pd.DataFrame, 
                             outfile=None, buff=3 ):
        
        """
        block_list : list of selected blocks' numbers\n
        G          : GridData object\n
        cc_df      : cutcell dataframe from cutcells_setup.dat\n
        outfile    : output file name. If None, using 'wall_var_projection.pkl'
        
        return     : dataframe of x-z plane data with coordinates and f\n
        
        Need data chunk with u,mu,wd ready.\n
        Only applicable to geometry with homogeneous shape along x
        """

# ----- initialize a list of dataframes containing groups of blocks' projection
        
        # a group of blocks defined as blocks sharing same x-z projection plane
        grouped_blocks = G.group_by_range( 'xz', block_list=block_list )
        f_visc_grp_ls = []   

# ----- loop over each group of blocks
     
        for group in grouped_blocks:
            
            # initialize a group of blocks' projection
            
            npx = G.g[group[0]-1].nx + buff*2
            npz = G.g[group[0]-1].nz + buff*2
            
            f_visc_grp = np.zeros( (npx,npz), dtype='f' )
            
# --------- loop over each block in the group
            
            for num in group:
            
# ------------- get cut cell dataframe and grid of this block 

                cc_df_block = cc_df[ cc_df['block_number'] == num ]
                cc_group = cc_df_block.groupby(['i','k'])
                
                g = G.g[num-1]
                
                npx = g.nx + buff*2
                npy = g.ny + buff*2
                npz = g.nz + buff*2
                
# ------------- prepare block data chunk

                # !!! i,j,k follow Fortran index, starting from 1.
                
                f_visc = np.zeros( (npx,npz), dtype='f' )
                
                data_df = self.bl[self.bl_nums.index(num)].df
                
                wd = np.array( data_df['wd'] ).reshape( npz, npy, npx )
                u  = np.array( data_df['u' ] ).reshape( npz, npy, npx )
                mu = np.array( data_df['mu'] ).reshape( npz, npy, npx )
                
# ------------- loop over each cut cell group
                
                for k in range( buff+1, g.nz+buff+1 ):
                    
                    # get a cut cell group dataframe
                    try:
                        df = cc_group.get_group((buff+1,k))
                    except:
                        # if there is no cut cell at this (x,z) location, skip and continue
                        # print(f"block {num}, point ({buff+1},{k}) no cut cell")
                        continue
                        
                    for i in range( buff+1, g.nx+buff+1 ):
                        
                        # when i == buff+1, find the interpolation stencil points
                        # (therefore, this code is just applicable to geometry
                        # with homogeneous shape in x direction.)
                        
                        if i == buff + 1:
                            
                            wd_cc  = [];   h         = []
                            y_cc   = [];   z_cc      = []
                            ny_cc  = [];   nz_cc     = []
                            fay    = [];   len_ratio = []
                            y_prj  = [];   z_prj     = []
                            jl_prj = [];   jr_prj    = []
                            kl_prj = [];   kr_prj    = []
                            
                            # loop over all cut cells sharing same x-z
                            
                            for cc_j in range( len(df) ):
                                
                                j = df['j'].iloc[cc_j]
                                
                                wd_cc.append( wd[k-1,j-1,i-1] )
                                y_cc .append( df['y'].iloc[cc_j])
                                z_cc .append( df['z'].iloc[cc_j])
                                ny_cc.append( df['ny'].iloc[cc_j])
                                nz_cc.append( df['nz'].iloc[cc_j])
                                
                                h.append( 0.50*np.sqrt((g.hy[buff]*ny_cc[cc_j])**2
                                                        +(g.hz[buff]*nz_cc[cc_j])**2))
                                
                                fay  .append( df['fay1'].iloc[cc_j]
                                            - df['fay0'].iloc[cc_j])
                                
                                len_ratio.append( ny_cc[cc_j] / np.sqrt( 
                                                    ny_cc[cc_j]**2 + nz_cc[cc_j]**2 ))
                                
                                # projection point coordinates
                                y_prj.append( y_cc[cc_j] 
                                            + (h[cc_j]-wd_cc[cc_j])*ny_cc[cc_j])
                                z_prj.append( z_cc[cc_j] 
                                            + (h[cc_j]-wd_cc[cc_j])*nz_cc[cc_j])
                                
                                # indices of interpolation stencil points
                                jl, jr = find_indices( g.gy, y_prj[cc_j] )
                                kl, kr = find_indices( g.gz, z_prj[cc_j] )
                                
                                jl_prj.append(jl)
                                jr_prj.append(jr)
                                kl_prj.append(kl)
                                kr_prj.append(kr)

                        # compute interpolated u, then friction of cut-cells 
                        # sharing the same x-z coordinate
                        
                        for cc_j in range( len(df) ):
                            
                            j = df['j'].iloc[cc_j]
                            
                            f = np.array([u[kl_prj[cc_j],jl_prj[cc_j],i],
                                          u[kr_prj[cc_j],jl_prj[cc_j],i],
                                          u[kl_prj[cc_j],jr_prj[cc_j],i],
                                          u[kr_prj[cc_j],jr_prj[cc_j],i]])
                            
                            u_prj = mth.bilin_interp( g.gz[kl_prj[cc_j]],
                                                      g.gz[kr_prj[cc_j]],
                                                      g.gy[jl_prj[cc_j]],
                                                      g.gy[jr_prj[cc_j]],
                                                      f,
                                                      z_prj[cc_j], y_prj[cc_j])                        
                            
                            # if geometry is fully 3D, extra term 4./3.u_norm 
                            # should be considered to add
                            
                            f_visc[i-1,k-1] += mu[k-1,j-1,i-1] * u_prj / h[cc_j] \
                                                * fay[cc_j] / len_ratio[cc_j]

# ------------- add frictions on each block of the group to compose the group's frictions

                print(f"block {num} has mean friction {np.mean(f_visc)}.\n")
                f_visc_grp += f_visc
                
# --------- match with coordinates
            
            xx,zz = np.meshgrid( g.gx, g.gz )
            
            df_fric_grp = pd.DataFrame(columns=['x','z','fric'])
            df_fric_grp['x'] = xx[buff:-buff,buff:-buff].flatten()
            df_fric_grp['z'] = zz[buff:-buff,buff:-buff].flatten()
            df_fric_grp['fric'] = f_visc_grp[buff:-buff,buff:-buff].T.flatten()
                
            f_visc_grp_ls.append( df_fric_grp )
            
# --------- save friction into StatisticsData.df_fric (a single pd.DataFrame)

        self.df_fric = pd.concat( f_visc_grp_ls )
        self.df_fric.reset_index( drop=True, inplace=True )
        self.df_fric.sort_values( by=['z','x'],inplace=True )
        
        if outfile is None: outfile = 'friction_projection.pkl'
        
        with open( outfile, 'wb' ) as f:
            pickle.dump( self.df_fric, f )


# ----------------------------------------------------------------------
# >>> wall variables projection                                    (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/26  - created
#
# Desc
#    - compute the pressure(fluctuation) distribution projected on x-z plane
#    - return a dataframe and output file
# ----------------------------------------------------------------------

    def wall_vars_projection( self,         block_list:list,
                              G:GridData,   cc_df:pd.DataFrame, 
                              outfile=None, buff=3 ):
        
        """
        block_list : list of selected blocks' numbers\n
        G          : GridData object\n
        cc_df      : cutcell dataframe from cutcells_setup.dat\n
        outfile    : output file name. If None, using 'wall_vars_projection.pkl'
        
        return     : dataframe of x-z plane data with coordinates and vars\n
        
        Need data chunk with <p> <pp> ready.\n
        Only applicable to geometry with homogeneous shape along x
        """

# ----- initialize a list of dataframes containing groups of blocks' projection
        
        grouped_blocks = G.group_by_range( 'xz', block_list=block_list )
        vars_grp_ls = []

# ----- loop over each group of blocks

        for group in grouped_blocks:
        
            # initialize a group of blocks' projection
            
            npx = G.g[group[0]-1].nx + buff*2
            npz = G.g[group[0]-1].nz + buff*2
            
            p_grp   = np.zeros( (npx,npz), dtype='f' )
            pp_grp  = np.zeros( (npx,npz), dtype='f' )
            mu_grp  = np.zeros( (npx,npz), dtype='f' )
            rho_grp = np.zeros( (npx,npz), dtype='f' )

# --------- loop over each block in the group
         
            for num in group:
                
# ------------- get cut cell dataframe and grid of this block 

                cc_df_block = cc_df[ cc_df['block_number'] == num ]
                cc_group = cc_df_block.groupby(['i','k'])
                
                g = G.g[num-1]
                
                npx = g.nx + buff*2
                npy = g.ny + buff*2
                npz = g.nz + buff*2
                
# ------------- prepare block data chunk

                # !!! i,j,k follow Fortran index, starting from 1.
                
                p_plane   = np.zeros( (npx,npz), dtype='f' )
                pp_plane  = np.zeros( (npx,npz), dtype='f' )
                mu_plane  = np.zeros( (npx,npz), dtype='f' )
                rho_plane = np.zeros( (npx,npz), dtype='f' )
                
                data_df = self.bl[self.bl_nums.index(num)].df
                
                p   = np.array( data_df['p' ]  ).reshape( npz, npy, npx )
                pp  = np.array( data_df['pp']  ).reshape( npz, npy, npx )
                mu  = np.array( data_df['mu']  ).reshape( npz, npy, npx )
                rho = np.array( data_df['rho'] ).reshape( npz, npy, npx )
                
# ------------- loop over each cut cell group

                for k in range( buff+1, g.nz+buff+1 ):
                    
                    # get a cut cell group dataframe
                    try:
                        df = cc_group.get_group((buff+1,k))
                    except:
                        # if there is no cut cell at this (x,z) location, skip and continue
                        #print(f"block {num}, point ({buff+1},{k}) no cut cell")
                        continue
                        
                    for i in range( buff+1, g.nx+buff+1 ):
                        
                        # when i == 1, find the interpolation stencil points
                        # (therefore, this code is just applicable to geometry
                        # with homogeneous shape in x direction.)
                        # Pressure and pressure fluctuation can use value on 
                        # cut cell directly, instead of interpolation.
                        
                        if i == buff + 1:
                            
                            fay    = []

                            # loop over all cut cells sharing same x-z
                            
                            for cc_j in range( len(df) ):
                                
                                fay.append( df['fay1'].iloc[cc_j] 
                                          - df['fay0'].iloc[cc_j])

                        # compute interpolated pressure and pressure fluctuation
                        
                        for cc_j in range( len(df) ):
                            
                            j = df['j'].iloc[cc_j]
                            
                            p_plane[i-1,k-1]   += p[k-1,j-1,i-1]   * fay[cc_j]
                            pp_plane[i-1,k-1]  += pp[k-1,j-1,i-1]  * fay[cc_j]
                            mu_plane[i-1,k-1]  += mu[k-1,j-1,i-1]  * fay[cc_j]
                            rho_plane[i-1,k-1] += rho[k-1,j-1,i-1] * fay[cc_j]

                print(f"block {num} has mean p {np.mean(p_plane)}")
                p_grp   += p_plane
                pp_grp  += pp_plane
                mu_grp  += mu_plane
                rho_grp += rho_plane

# --------- match with coordinates
            
            xx,zz = np.meshgrid( g.gx, g.gz )
            
            df_wall_grp = pd.DataFrame(columns=['x','z','p','pp','mu','rho'])
            df_wall_grp['x']   =        xx[buff:-buff,buff:-buff].flatten()
            df_wall_grp['z']   =        zz[buff:-buff,buff:-buff].flatten()
            df_wall_grp['p']   =     p_grp[buff:-buff,buff:-buff].T.flatten()
            df_wall_grp['pp']  =    pp_grp[buff:-buff,buff:-buff].T.flatten()
            df_wall_grp['mu']  =    mu_grp[buff:-buff,buff:-buff].T.flatten()
            df_wall_grp['rho'] =   rho_grp[buff:-buff,buff:-buff].T.flatten()
            
            # compute pressure fluctuation
            df_wall_grp['p`'] = np.sqrt( np.array( df_wall_grp['pp']-df_wall_grp['p']**2 ))
            
            vars_grp_ls.append( df_wall_grp )
            
# --------- save wall vars into StatisticsData.df_wall (a single pd.DataFrame)

        self.df_wall = pd.concat( vars_grp_ls )
        self.df_wall.reset_index( drop=True, inplace=True )
        self.df_wall.sort_values( by=['z','x'],inplace=True )
        
        if outfile is None: outfile = 'wall_vars_projection.pkl'
        
        with open( outfile, 'wb' ) as f:
            pickle.dump( self.df_wall, f )


# ----------------------------------------------------------------------
# >>> extract wall varibles of smooth wall                     (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def extract_wall_vars_sw( self,       block_list:list,
                              G:GridData, outfile=None,
                              buff=3):
        
        """
        block_list : list of selected blocks' numbers\n
        G          : GridData object\n  
        outfile    : output file name. If None, using 'wall_vars_projection.pkl'      
        
        return     : dataframe of x-z plane data with coordinates and vars\n
        vars      : ['x','z','p','pp','mu','rho','fric','p`']\n
        
        Only applied to the smooth wall case. Extracting wall variables.
        Need data chunk with <p>,<pp>,<u>,<mu>,<rho> ready.
        Only applicable to smooth wall case.
        """
        
        for num in block_list:

# --------- get grid of this block

            g = G.g[num-1]
            
            npx = g.nx + buff*2
            npy = g.ny + buff*2
            npz = g.nz + buff*2
            
# --------- prepare block data chunk

            data_df = self.bl[self.bl_nums.index(num)].df
            
            p   = np.array( data_df['p' ]  ).reshape( npz, npy, npx )
            pp  = np.array( data_df['pp']  ).reshape( npz, npy, npx )
            mu  = np.array( data_df['mu']  ).reshape( npz, npy, npx )
            rho = np.array( data_df['rho'] ).reshape( npz, npy, npx )
            u   = np.array( data_df['u']   ).reshape( npz, npy, npx )
            
# --------- get slice from those data chunks
            
            dy = g.gy[buff]
            
            p_plane   = p[:,buff,:]
            pp_plane  = pp[:,buff,:]
            mu_plane  = mu[:,buff,:]
            rho_plane = rho[:,buff,:]
            u_plane   = u[:,buff,:]
            f_visc    = u_plane * mu_plane / dy

            xx,zz = np.meshgrid( g.gx, g.gz)
            
            df_wall = pd.DataFrame(columns=['x','z','p','pp','fric'])
            df_wall['x']    =       xx[buff:-buff,buff:-buff].flatten()
            df_wall['z']    =       zz[buff:-buff,buff:-buff].flatten()
            df_wall['p']    =  p_plane[buff:-buff,buff:-buff].flatten()
            df_wall['pp']   = pp_plane[buff:-buff,buff:-buff].flatten()
            df_wall['mu']   =  mu_plane[buff:-buff,buff:-buff].flatten()
            df_wall['rho']  = rho_plane[buff:-buff,buff:-buff].flatten()
            df_wall['fric'] =   f_visc[buff:-buff,buff:-buff].flatten()
            
            # compute pressure fluctuation
            df_wall['p`'] = np.sqrt( np.array( df_wall['pp']-df_wall['p']**2 ))

            self.bl[self.bl_nums.index(num)].df_wall = df_wall

            p_ave      = np.mean( np.array(df_wall['p'])  )
            p_fluc_ave = np.mean( np.array(df_wall['p`']) )
            fric_ave   = np.mean( np.array(df_wall['fric']) )
            
            print(f"block {num} has mean p {p_ave},",end='')
            print(f" p` {p_fluc_ave}, friction {fric_ave}.\n")

# --------- save wall vars into StatisticsData.df_wall (a single pd.DataFrame)

        self.df_wall = pd.concat([self.bl[self.bl_nums.index(num)].df_wall for num in block_list])
        self.df_wall.reset_index( drop=True, inplace=True )
        self.df_wall.sort_values( by=['z','x'], inplace=True )
        
        if outfile is None: outfile = 'wall_vars_projection.pkl'
        
        with open( outfile, 'wb' ) as f:
            pickle.dump( self.df_wall, f )
            

# ----------------------------------------------------------------------
# >>> assign_wall_dist                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/04/25  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def assign_wall_dist( self, block_list, wd_snap ):
        
        """
        block_list : list of selected blocks' numbers
    
        wd_snap : Snapshot instance with wall distance field
        
        copy wall distance field from wd_snap to self.bl[].df
        """
        
        for num in block_list:
            
            if self.bl[self.bl_nums.index(num)].num != wd_snap.snap_data[wd_snap.bl_nums.index(num)].num:
                raise ValueError("block number not match when assigning walldist.")
            
            self.bl[self.bl_nums.index(num)].df['wd'] = wd_snap.snap_data[wd_snap.bl_nums.index(num)].df['wd']


# ----------------------------------------------------------------------
# >>> compute bubble volume                                         (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/04/25  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_bubble_volume( self, G:GridData, block_list,
                               cc_df=None, roughwall=False, buff=3,
                               y_threshold=None ):
        
        """
        G     : GridData instance
        cc_df : cutcell dataframe from cutcells_setup.dat
        
        return: separation bubble volume
        Need data chunk with u ready.
        G should contain cell volume.
        wd (wall distance) should be contained in self.bl[num-1].df
        y_threshold : only calculate bubble volume above this y threshold
        """
        
        vol_bubble     = 0.0
        vol_bubble_thr = 0.0
        
        for num in block_list:
            
            df = self.bl[self.bl_nums.index(num)].df
            g  = G.g[num-1]

            npx = g.nx + buff*2
            npy = g.ny + buff*2
            npz = g.nz + buff*2
            
            u = np.array( df['u'] ).reshape( npz, npy, npx )
            
            if y_threshold is not None:
                y  = np.meshgrid(g.gx,g.gy,g.gz, indexing='ij')[1]
                identifier_thr = (u<0.0) & (y.T>y_threshold)
                identifier_thr = identifier_thr*1.0
                
            identifier = u < 0.0
            identifier = identifier*1.0
            
            vol = g.vol
            
            if roughwall:
                
                temp_df   = cc_df[ cc_df['block_number'] == num ]
                wall_dist = np.array( df['wd'] )
                g.assign_vol_fra( df=temp_df, wall_dist=wall_dist )
            
            else:
                
                g.assign_vol_fra()
            
            vol_bubble_block = vol*identifier*(g.vol_fra.T)
            vol_bubble      += np.sum(vol_bubble_block[buff:-buff,buff:-buff,buff:-buff])
            
            if y_threshold is not None:
                vol_bubble_block_thr = vol*identifier_thr*(g.vol_fra.T)
                vol_bubble_thr      += np.sum(vol_bubble_block_thr[buff:-buff,buff:-buff,buff:-buff])
            
        self.vol_bubble = vol_bubble
        
        if y_threshold is not None:
            self.vol_bubble_thr = vol_bubble_thr
            return vol_bubble, vol_bubble_thr
        
        else:
            return vol_bubble


# ----------------------------------------------------------------------
# >>> spanwise average                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/03/14  - created
#
# Desc
#     Do spanwise average and store the averaged results into a pickeled 
#     StatisticData (inca's statistic.bin requires all dozens of variables, 
#     too large and unnecessary).
# 
# ----------------------------------------------------------------------

    def spanwise_average( self, blocklist, vars, 
                          buff=3, rescale=[0.0,0.0,0.0,1.0,1.0,1.0] ):
        
        """
        blocklist : blocks that apply spanwise average
        vars      : list of variables
        """
        
        grid3d    = self.grid3d
        stat_type = self.stat_type
        
        if stat_type != 'block':
            raise ValueError("Only 'block' type is supported for spanwise average.")
        
        # - drop ghost cells
        
        self.drop_ghost( blocklist, buff=3, mode='symmetry' )
        
        # - get list of grouped block
        grouped_block_list = grid3d.group_by_range('xy', block_list=blocklist)
        
        vtk_blocks = list()
        
        # -- loop over the list
        for group in grouped_block_list:
            
            g  = grid3d.g[group[0]-1]
            px = (g.px[buff:-buff]+rescale[0])
            py = (g.py[buff:-buff]+rescale[1])
            pz = np.array([-0.1,0.1])
            
            bl_vtk = create_3d_vtkRectilinearGrid( px/rescale[3], py/rescale[4], pz/rescale[5] )
            
            for var in vars:
                
                var_data = list()
                
                for num in group:
                    bl_data = self.bl_clean[self.bl_nums_clean.index(num)].df[var]
                    var_data.append( np.array(bl_data).reshape(g.nz,g.ny,g.nx).mean(axis=0) )
                
                var_data = np.array(var_data).mean(axis=0)
                
                if   var_data.size != (len(px)-1)*(len(py)-1):
                    raise ValueError("Data length does not match.")
                elif var_data.size == 0:        
                    raise ValueError(f"There is no data in blocks {group}.")

                bl_vtk = add_var_vtkRectilinearGrid( bl_vtk, var, var_data )

            vtk_blocks.append( bl_vtk )
        
        dataset = create_multiblock_dataset( vtk_blocks )
        
        return dataset


# ----------------------------------------------------------------------
# >>> spanwise periodic average                                 (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/02/25  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def spanwise_periodic_average( self, blocklist, vars, period ):
        
        """
        Do periodic average and store the averaged results to original bl.df
        
        blocklist : blocks that apply periodic average
        vars      : list of variables
        period    : length of period
        """
        
        grid3d    = self.grid3d
        stat_type = self.stat_type
        
        # - get list of grouped block
        grouped_block_list = grid3d.group_by_range('xy', block_list=blocklist)
        
        # - check if x,y,z is added to the dataframe of blocks
        z_is_ready = 'z' in self.bl[self.bl_nums.index(blocklist[0])].df.columns
        x_is_ready = 'x' in self.bl[self.bl_nums.index(blocklist[0])].df.columns
        
        if not (x_is_ready or z_is_ready):
            raise ValueError("Please match grid first and add coordinates to the dataframes of blocks.")
        
        # - loop over the list
        for group in grouped_block_list:
        
        # ---- check if domain size is multiple of period
            z0 = min([grid3d.g[num-1].lz0 for num in group])
            z1 = max([grid3d.g[num-1].lz1 for num in group])
            if np.round((z1-z0) % period,7) != 0.0:
                raise ValueError(f"Domain size in z is not a multiple of period.")
       
        # ---- add coordiante remainder to the dataframe 
        
            for num in group:
                df = self.bl[self.bl_nums.index(num)].df
                
                # compute remainder of z
                df['remainder'] = np.round(np.array(df['z']) % period,7)
            
            # - concat all dataframe from blocks
            group_df = pd.concat([self.bl[self.bl_nums.index(num)].df for num in group ])
            
            for var in vars:
                if stat_type == 'block':
                    group_df[var] = group_df.groupby(['remainder','x','y'])[var].transform('mean')
                elif stat_type == 'Y':
                    group_df[var] = group_df.groupby(['remainder','x'])[var].transform('mean')
                elif stat_type == 'X':
                    group_df[var] = group_df.groupby(['remainder','y'])[var].transform('mean')
            
            # - distribute averaged data back to each dataframe

            k = 0
            for num in group:
                df     = self.bl[self.bl_nums.index(num)].df
                for var in vars:
                    df[var] = np.array(group_df[var][k:k+len(df)])
                k += len(df)


# ----------------------------------------------------------------------
# >>> create vtk multiblock dataset (statistics)      (Nr.)
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

    def create_vtk_multiblock( self, block_list, vars, buff=3, 
                              mode='symmetry', rescale=[0.0,0.0,0.0,1.0,1.0,1.0] ):
        
        """
        write statistics into vtm file (multiblock vtk)
        
        block_list : list of selected blocks' numbers
        vars       : list of variables to be written
        buff       : number of ghost layers\n
        mode       : 'symmetry' or 'oneside'\n
        rescale    : rescale the grid points. [x_shift, y_shift, z_shift, x_norm, y_norm, z_norm]
        """

# ----- check if grid data is ready

        if self.grid3d is None:
            raise ValueError("Please set self.grid3d before read_stat")


# ----- drop ghost cells

        self.drop_ghost( block_list, buff=buff, mode=mode )
        
        if mode == 'symmetry':  buffl = buff; buffr = buff
        elif mode == 'oneside': buffl = buff; buffr = buff-1

# ----- setup vtk file

        vtk_blocks =  list()
        
        for bl in self.bl_clean:
            
            # only write selected blocks, different to SnapData
            
            if bl.num not in block_list:
                continue
            
            bl_num = bl.num
            g = bl.g
            
            px = (g.px[buffl:-buffr]+rescale[0])/rescale[3]
            py = (g.py[buffl:-buffr]+rescale[1])/rescale[4]
            pz = (g.pz[buffl:-buffr]+rescale[2])/rescale[5]
                        
            # build one vtk block
            
            if bl.npx == 1: px = np.array([0.0])
            if bl.npy == 1: py = np.array([0.0])
            if bl.npz == 1: pz = np.array([0.0])
            
            bl_vtk = create_3d_vtkRectilinearGrid( px, py, pz )
            
            for var in vars:
                
                var_data = np.array(bl.df[var])
                
                if len(var_data) != bl.npx*bl.npy*bl.npz:
                    raise ValueError(f"Data length not match for variable {var} in block {bl_num}.")
                elif len(var_data) == 0:
                    raise ValueError(f"Data length is zero for variable {var} in block {bl_num}.")

                bl_vtk = add_var_vtkRectilinearGrid( bl_vtk, var, var_data )
                
            vtk_blocks.append( bl_vtk )
        
        # build the multiple blocks dataset     
        dataset = create_multiblock_dataset(vtk_blocks)
        
        return dataset 


# ----------------------------------------------------------------------
# >>> write statistics into vtm (multiblock vtk) file        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/12  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def write_vtm( self, filename, vars, block_list, buff=3 ):
        
        """
        write statistics into vtm file (multiblock vtk)
        
        filename   : filename of output statistics
        vars       : list of variables to be written
        block_list : list of selected blocks' numbers
        """
        
        # build the multiple blocks dataset     
        dataset = self.create_vtk_multiblock( block_list, vars, buff=buff )
        
        write_vtm_file( filename, dataset )
                    

# ----------------------------------------------------------------------
# >>> get a slice from 3D Statistic                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/25  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def get_slice( self, slic_type, loc:float, buff=3 ):
        
        """
        slic_type : 'X','Y', or 'Z'
        loc       : location of slice
        
        return    : 2D statistic instance
        """

# ----- check if current statistic is 3D?

        tbl = self.bl[0]
        if tbl.npx==1 or tbl.npy==1 or tbl.npz==1:
            raise TypeError("Current statistic is not 3D data.")

# ----- check if the grid file is available

        if self.grid3d is None:
            raise ValueError("Please read in grid file first!")
        
        else: grid3d = self.grid3d   


# ----- select the blocks that intersect the plane
#       1. intersected blocks' number 2. indexes of sliced location

        bl_intersect = []
        indx_slic    = []

        # loop over all blocks within this 3d statistic
        
        for block in self.bl:
            
            bl_num = block.num
            grd    = grid3d.g[bl_num-1]
            
            if slic_type == 'X':
                
                # if intersect, keep this block number and 
                # find the index of slice location
                if loc >= grd.lx0 and loc < grd.lx1:
                    
                    bl_intersect.append( bl_num )
                
                    for i in range( buff, grd.nx+buff ):
                        if (grd.gx[i]-0.5*grd.hx[i] <= loc < grd.gx[i]+0.5*grd.hx[i]):
                            
                            indx_slic.append( i )
                            break
            
            elif slic_type == 'Y':
                
                if loc >= grd.ly0 and loc < grd.ly1:
                    
                    bl_intersect.append( bl_num )
                    
                    for j in range( buff, grd.ny+buff ):
                        if (grd.gy[j]-0.5*grd.hy[j] <= loc < grd.gy[j]+0.5*grd.hy[j]):
                            
                            indx_slic.append( j )
                            break
            
            elif slic_type == 'Z':
                
                if loc >= grd.lz0 and loc < grd.lz1:
                    
                    bl_intersect.append( bl_num )
                    
                    for k in range( buff, grd.nz+buff ):
                        if (grd.gz[k]-0.5*grd.hz[k] <= loc < grd.gz[k]+0.5*grd.hz[k]):
                            
                            indx_slic.append( k )
                            break
        
# ----- check if any block intersected

        if len(bl_intersect) == 0:
            raise ValueError("No block is sliced! Check slice location.")
        
        if self.verbose:
            print(f" {len(bl_intersect)} blocks are sliced.\n")
            print(" index of block  --  index of cut cell")
            for i in range( len(bl_intersect) ):
                print(f" {bl_intersect[i]:15d} -- {indx_slic[i]:12d}")

# ----- init a slice statistic instance

        stat_2d = StatisticData()
        
        # get sliced data from 3d statistic and fill in a 2d statistic
        
        for block in self.bl:
            
            bl_num = block.num
            
            if bl_num in bl_intersect:
                
                grd = grid3d.g[bl_num-1]
                idx = indx_slic[ bl_intersect.index(bl_num) ]
                
                N1 = grd.nx + buff*2
                N2 = grd.ny + buff*2
                N3 = grd.nz + buff*2
                
                vars   = block.df.columns
                n_var  = len(vars)
                df_sol = deepcopy( block.df.values)
                sol    = np.array(df_sol).T.reshape(n_var,N3,N2,N1)
                
                if slic_type == 'X':
                    
                    dims = [1,N2,N3]
                    sol  = sol[:,:,:,idx].reshape(n_var,N3*N2)
                
                elif slic_type == 'Y':
                    
                    dims = [N1,1,N3]
                    sol  = sol[:,:,idx,:].reshape(n_var,N3*N1)
                
                elif slic_type == 'Z':
                    
                    dims = [N1,N2,1]
                    sol  = sol[:,idx,:,:].reshape(n_var,N2*N1)
                    
                df_sol   = pd.DataFrame(sol.T, columns=vars)
                
                bl_slice = BlockData()
                bl_slice.fill_with_data( bl_num, dims, df_sol )
                stat_2d.bl.append( bl_slice )
                
# ----- fill in headers

        if self.format == 1:
            
            headers = ['format',         'n_samples',        'sample_step',  
                       'start_step',     'sample_time',      'start_time', 
                       'meanvalues',     'doublecorr',       'triplecorr',
                       'quadruplecorr',  'autocorr',         'mean_invar',
                       'schlieren',      'cavitation_stats', 'vapor_gas_stats',
                       'rste',           'thermo',           'visc_diff'
                       ]
        
        if self.format == 2:

            headers = ['format',         'n_samples',        'sample_step',  
                       'start_step',     'npv',              'nscalars',
                       'nspecies',       'sample_time',      'start_time', 
                       'meanvalues',     'doublecorr',       'triplecorr',
                       'quadruplecorr',  'autocorr',         'mean_invar',
                       'schlieren',      'cavitation_stats', 'vapor_gas_stats',
                       'rste',           'thermo',           'visc_diff'
                       ]

        for attr, value in self.__dict__.items():
            if attr in headers:
                setattr( stat_2d, attr, deepcopy(value) )

# ----- return 2d statistic instance

        return stat_2d



# ----------------------------------------------------------------------
# >>> write statistics                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/28  - created
#
# Desc
#     write out statistics.bin from a StatisticData instance
# ----------------------------------------------------------------------

    def write_statistic( self, filename ):

        """
        filename: output filename
        """
        
        verbose = self.verbose
        
# ----- write statistics header

        # default size
        sin = 4; slg = 4; sfl = 8
        
        with open( filename, 'wb' ) as f:
            
            # header format
            
            write_int_bin( self.format,      f, sin )
            
            # number of samples, sample step & start step
            
            write_int_bin( self.n_samples,   f, sin )
            write_int_bin( self.sample_step, f, sin )
            write_int_bin( self.start_step,  f, sin )
            
            # check format
            
            if self.format == 2:
                
                # number of transported variables
                
                write_int_bin( self.npv,      f, sin )
                write_int_bin( self.nscalars, f, sin )
                write_int_bin( self.nspecies, f, sin )
        
            # sample time & start time
            
            write_flt_bin( self.sample_time, f, sfl )
            write_flt_bin( self.start_time,  f, sfl )
            
            buf_log = [ self.meanvalues,
                        self.doublecorr,
                        self.triplecorr,
                        self.quadruplecorr,
                        self.autocorr,
                        self.mean_invar,
                        self.schlieren,
                        self.cavitation_stats,
                        self.vapor_gas_stats,
                        self.rste,
                        self.thermo,
                        self.visc_diff,
                        ]
            
            write_log_bin( buf_log, f, slg )

            if verbose: print("Finish writing statistics header.")
        
# --------- write statistics body
            
            for bl in self.bl:
                
                # block number
                
                write_int_bin( bl.num, f, sin )

                # dimensions of grids, N1, N2, N3
                
                write_int_bin( bl.npx, f, sin )
                write_int_bin( bl.npy, f, sin )
                write_int_bin( bl.npz, f, sin )
                
                # sol(n_var, N3, N2, N1)
                
                write_flt_bin( np.array(bl.df.values).T, f, sfl )

        print(f"Finish writing {filename}.")

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/28  - created
#
# Desc
#
# ----------------------------------------------------------------------
    
    @property
    def full_vars( self ):

        mean_ls = ['u','v','w','rho','rhoE','p','T','mu']
        
        cor2_ls = ['uu','uv','uw','urho'  ,'urhoE'   ,'up'   ,'uT'   ,'umu'   ,
                        'vv','vw','vrho'  ,'vrhoE'   ,'vp'   ,'vT'   ,'vmu'   ,
                             'ww','wrho'  ,'wrhoE'   ,'wp'   ,'wT'   ,'wmu'   ,
                                  'rhorho','rhorhoE' ,'rhop' ,'rhoT' ,'rhomu' ,
                                           'rhoErhoE','rhoEp','rhoET','rhoEmu',
                                                      'pp'   ,'pT'   ,'pmu'   ,
                                                              'TT'   ,'Tmu'   ,
                                                                      'mumu'  ]
        
        return mean_ls + cor2_ls

# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/11/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

# -------- test get_slice ----------------

    path  = '/home/wencanwu/my_simulation/temp/220927_lowRe/results'
    gfile = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/inca_grid.bin'
    
    os.chdir( path )
    
    grid = GridData( gfile )
    grid.read_grid()
    
    blocklist, index = grid.select_sliced_blockgrids( 'Z', 0.1 )
    
    print( blocklist )
    
    stat3d = StatisticData( path+'/statistics.bin' )
    stat3d.verbose = True
    stat3d.read_statistic( block_list=blocklist, vars_in = stat3d.full_vars )
    stat3d.grid3d = grid
    
    stat2d = stat3d.get_slice( 'Z', 0.1 )
    stat2d.verbose =True 
    
    stat2d.write_statistic( 'slicez.bin' )
    
    stattest = StatisticData( 'slicez.bin' )
    
    stattest.verbose = True
    
    stattest.read_statistic( block_list=blocklist, vars_in=stattest.full_vars )


# -------- test write vtm ----------------

    # path = "/media/wencan/Expansion/temp/231124/results"
    # outpath = "/media/wencan/Expansion/temp/231124/vtk"
    
    # os.chdir(path)
    
    # box = [-999,999,-999,999,-999,999]
    
    
    # G = GridData( 'inca_grid.bin')
    # G.read_grid()
    # block_list =  G.select_blockgrids( box )
    
    
    # S = StatisticData( 'statistics.bin' )
    
    # with open( 'statistics.bin', 'rb') as f:
        
    #     S.verbose = True
        
    #     vars = ['u','v','w','p','rho']
        
    #     S.read_stat_header( f )
    #     S.read_stat_body( f, block_list, vars)
    
        
    # S.match_grid( block_list, G, add_to_df=False)
    # S.grid3d = G
    
    # os.chdir(outpath)
    # S.write_vtm( 'statistics.vtm', vars, block_list )
