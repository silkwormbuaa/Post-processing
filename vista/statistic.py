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

import numpy             as     np

import pandas            as     pd

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

    def __init__( self, stat_file ):
        
        # directory to the statistic file
        
        self.file = stat_file

        # file size
        
        self.fsize = os.stat(self.file).st_size
        
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

    def read_stat_body( self, file , fill, vars ):
        
        end_of_file = False
               
        # New position after reading headers
        
        file.seek( self.header_size )
        
        self.pos = self.header_size
        
        # vars determine the index of data sequences to read
        
        vars_indx = self.vars_to_indx( vars )
        
        while not end_of_file:
        
            # read in block one by one
            # only blocks in fill list will be filled with data chunk
        
            self.bl.append( BlockData(file, self.n_var, fill, vars, vars_indx) )
            
            self.pos = self.pos + self.bl[-1].size
                                    
            if self.pos >= self.fsize: end_of_file = True


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
        
        # list of variables for inquiry
        
        mean_ls = ['u','v','w','rho','rhoE','p','T','mu']
        
        cor2_ls = ['uu','uv','uw','urho'  ,'urhoE'   ,'up'   ,'uT'   ,'umu'   ,
                        'vv','vw','vrho'  ,'vrhoE'   ,'vp'   ,'vT'   ,'vmu'   ,
                             'ww','wrho'  ,'wrhoE'   ,'wp'   ,'wT'   ,'wmu'   ,
                                  'rhorho','rhorhoE' ,'rhop' ,'rhoT' ,'rhomu' ,
                                           'rhoErhoE','rhoEp','rhoET','rhoEmu',
                                                      'pp'   ,'pT'   ,'pmu'   ,
                                                              'TT'   ,'mumu'  ]   
        
        # index list of vars
        
        indx = []
        
        # from variables list to corresponding indexes list
        
        for var in vars:
            
            displace = 0
            
            if var in mean_ls:
                
                indx.append( displace + mean_ls.index(var) )

            if self.meanvalues: displace += 8
            
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
    
    def match_grid( self, G, block_list ):
        
        for num in block_list:
            
            g = G.g[num-1]
            
            X,Y,Z = np.meshgrid( g.gx, g.gy, g.gz, indexing='ij' )
            
            self.bl[num-1].df['x'] = X.T.flatten()
            self.bl[num-1].df['y'] = Y.T.flatten()
            self.bl[num-1].df['z'] = Z.T.flatten()
            
            # adding vol_fra !! original vol_fra has i,j,k order, 
            # should be transpose as k,j,i
            
            self.bl[num-1].df['vol_frac'] = np.ravel( G.g[num-1].vol_fra.T )
            

# ----------------------------------------------------------------------
# >>> drop ghost cells                                           (Nr.)
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

    def drop_ghost( self, G, block_list, buff=3 ):
        
        """
        Applicable to 3D statistics data.
        G : GridData instance
        block_list : list of blocks that are going to drop ghost cells
        
        return : self.bl[num-1].df
        """
        
        for num in block_list:
            
            npx = G.g[num-1].nx + buff*2
            npy = G.g[num-1].ny + buff*2
            npz = G.g[num-1].nz + buff*2
            
            vars = self.bl[num-1].df.columns
            
            data_chunk = None
            
            for var in vars:
                
                data = np.array( self.bl[num-1].df[var] )
                data = data.reshape( npz, npy, npx )
                data = data[buff:npz-buff,buff:npy-buff,buff:npx-buff].flatten()
                
                if data_chunk is None: data_chunk = [data]
                else: data_chunk.append(data)
            
            data_chunk = np.array(data_chunk).T
            
            self.bl[num-1].df = pd.DataFrame(data_chunk,columns=vars)
            
        #    print( self.bl[num-1].df)    
                


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
        block_list: list of blocks that are going to compute new variables
        vars_new: list of str representing new vars 
                  ['mach','grad_rho','grad_RS']
        """
        
        
        for num in block_list:
            
            df = self.bl[num-1].df
            
# --------- compute Mach number 

            if "mach" in vars_new:
                
                u = np.array( df['u'] )
                v = np.array( df['v'] )
                w = np.array( df['w'] )
                T = np.array( df['T'] )
                gamma = 1.40
                R = 287.0508571
                
                mach = np.sqrt( u*u+v*v+w*w ) / np.sqrt( gamma*R*T )
                
                self.bl[num-1].df['mach'] = mach
            
# ---------- compute density gradient

            if "grad_rho" in vars_new:
                
                pass
            
# ---------- compute Reynolds Stress gradients

            if "grad_RS" in vars_new:
                
                pass
                
                
                
                
        


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

    def compute_profile( self, block_list, bbox, vars, RS=True, outfile=False ):
        
        """
        block list: selected block list (within bounding box)
        bbox: bounding box coordinates [xmin,xmax,ymin,ymax,zmin,zmax]
        vars: variable names list
        RS: set True to compute Reynolds Stresses, otherwise drop them
        outfile: assign outfile name or use default one ''
        """
        
        
# ----- collect data frame from all filled blocks
        
        df = pd.concat( [self.bl[num-1].df for num in block_list] )
        

        # reset indexes in case repeated indexes from different blocks
        
        df.reset_index( drop=True, inplace=True )
        
        print(df)
        
        
# ----- slim down data with bbox
        
        df.drop( df[ (df['x'] < bbox[0]) | (df['x'] > bbox[1]) |
                     (df['y'] < bbox[2]) | (df['y'] > bbox[3]) |
                     (df['z'] < bbox[4]) | (df['z'] > bbox[5]) ].index,
                 inplace=True )
        
        print(df)
        
        # compute turbulence (Reynolds Stress)
        
        if RS:
            df['uu'] = np.array(df['uu']) - np.array(df['u'])*np.array(df['u'])
            df['vv'] = np.array(df['vv']) - np.array(df['v'])*np.array(df['v'])
            df['ww'] = np.array(df['ww']) - np.array(df['w'])*np.array(df['w'])
            df['uv'] = np.array(df['uv']) - np.array(df['u'])*np.array(df['v'])
            
            print(df)
            
# ----- do averaging in x-z plane
        
        ys = np.sort( np.unique( np.array( df['y'] )))
        
        data_chunk = None
        
        for y in ys:
            
            df_temp = df[ df['y']==y ]
            
            vol_frac = np.array( df_temp['vol_frac'] )
            vol_total = np.sum( vol_frac )
            if vol_total < 0.0000001: vol_total = float('inf')
            
            buff = [y]
            
            for var in vars:
                
                v = np.sum( np.array(df_temp[var])*vol_frac ) / vol_total
                buff.append(v)
            
                
            if data_chunk is None: data_chunk = [buff]
            else: data_chunk.append( buff )
            
        
        data_chunk = np.array(data_chunk).reshape( len(ys), len(vars)+1 )
        
        df_profile = pd.DataFrame(data_chunk,columns=['y']+vars)
        
        pd.set_option('display.max_rows', None)  # 显示所有行
        
        print( df_profile )
        
# ----- drop points below wall and set values(except p,rho,T) at wall manually
        
        df_profile.drop( df_profile[ df_profile['p']==0.0 ].index, inplace=True)
        
        df_profile.reset_index( drop=True, inplace=True )
        
        for var in vars:
            if var not in ['p','rho','T']:
                df_profile.loc[0,var] = 0.0
        
# ----- output profile into txt

        if outfile is False: outfile = 'profile_spanwisemean.dat'

        df_profile.to_string(outfile, index=False,
                            float_format='%15.7f',
                            justify='left')
         
            

# ----------------------------------------------------------------------
# >>> get slice                                                (Nr.)
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

    def get_slice( self, G, block_list, indx_slic, slic_type, buff=3 ):
        
        """
        G : corresponding GridData instance
        block_list: list of sliced blocks
        indx_slic: list of slice indexes on each block
        slic_type : 'X','Y' or 'Z'
        loc : location coordinate
        
        return : dataframe of sliced results
        """
        
# ----- should somehow check if Statistic is 3D or 2D data?
        
        # pass
        
        for i, bl_num in enumerate( block_list ):

# ----- do slicing on each bl.df
            
            bl   = self.bl[bl_num-1]
            g    = G.g[bl_num-1]
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

# ----- compute gradients if necessary

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

        df_slice = pd.concat( [self.bl[num-1].df for num in block_list] )
        
        df_slice.reset_index( drop=True, inplace=True)
        
        return df_slice