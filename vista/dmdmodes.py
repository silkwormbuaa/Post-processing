# -*- coding: utf-8 -*-
'''
@File    :   dmdmodes.py
@Time    :   2023/06/05 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import pickle
import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import griddata

from   .plot_style       import plot_dmd_mode 
from   .io_vtk           import create_multiblock_dataset
from   .io_vtk           import add_var_vtkRectilinearGrid
from   .io_vtk           import create_3d_vtkRectilinearGrid

class DMDMode:
# ----------------------------------------------------------------------
# >>> Class of DMD mode                                           (Nr.)
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
# ----------------------------------------------------------------------

    def __init__( self, modefile = None ):
        
        '''
        indx:       index of the mode among all standard DMD modes
        alpha:      amplitude from standard DMD
        alpha_sp:   amplitude before polishing from SPDMD 
        alpha_pol:  amplitude after polishing from SPDMD
        St:         Strouhal number of modes
        mu:         eigenvalues from standard DMD
        '''
        
        self.indx      = None
        self.alpha     = None
        self.alpha_sp  = None
        self.alpha_pol = None
        self.mu        = None
        self.St        = None
        self.Phi       = None

        if modefile is not None:
            
            self.read_mode_file( modefile )
        


# ----------------------------------------------------------------------
# >>> Read mode file                                             (Nr.)
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
# ----------------------------------------------------------------------

    def read_mode_file( self, filename ):
        
        """
        filename: dmd_modes binary file
        """
        
        with open( filename, 'rb' ) as f:
            
            self.indx      = pickle.load( f )
            
            self.alpha     = pickle.load( f )
            
            self.alpha_sp  = pickle.load( f )
            
            self.alpha_pol = pickle.load( f )
            
            self.mu        = pickle.load( f )
            
            self.St        = pickle.load( f )
            
            self.Phi       = pickle.load( f )
        


class DMDModes:
# ----------------------------------------------------------------------
# >>> a group of DMD modes                                        (Nr.)
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
# ----------------------------------------------------------------------

    def __init__( self ):
        
        """
        Initialize a void DMDModes
        """
        
        self.modes = []
        
        self.grids_interp = None
        
        self.case_parameters = None

        self.recons_data = None
        
        self.recons_std_dmd = None


# ----------------------------------------------------------------------
# >>> Add mode to group                                          (Nr.)
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
# ----------------------------------------------------------------------

    def add_mode( self, mode ):
        
        """
        mode: an instance of DMDMode()
        """
        
        self.modes.append( mode )
        


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
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
# ----------------------------------------------------------------------

    def reconstruct( self, step ):
        
        """
        step: advance how many dt. \n 
        Set to 1 to reconstruct initial flow field. 
        """
        
        self.step = step
        
        self.indxes = [mode.indx for mode in self.modes]

        self.Phis = np.array( [mode.Phi for mode in self.modes] ).T
        
        self.alphas = [mode.alpha for mode in self.modes]

        self.alphas_pol = [mode.alpha_pol for mode in self.modes]
        
        self.mus = [mode.mu for mode in self.modes]
        
        self.Sts = [mode.St for mode in self.modes] 

        # build a new dataframe of above variables
        
        df_ind = pd.DataFrame( self.indxes, columns=['indxes'])
        
        df_ind['Sts'] = np.array( self.Sts )
        df_ind['alphas'] = np.array( self.alphas )
        df_ind['alphas_pol'] = np.array( self.alphas_pol )
        df_ind['mus'] = np.array(self.mus)
        
        df_ind = df_ind.drop(df_ind[df_ind['Sts'] < 0].index)
        df_ind = df_ind.sort_values(by='Sts')

        self.df_ind = df_ind
        
        # reconstruct and time advance
        
        self.vand = np.vander( self.mus, self.step, increasing=True )

        if not self.modes:
        
            raise ValueError("No mode is available!")
        
        else:
            
            self.recons_data = self.Phis @ np.diag(self.alphas_pol) @ self.vand

        # clear the list of modes
        
        self.modes.clear()



# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/12/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def create_vtk_multiblock( self, vars, block_list, snap_type, grid3d,
                               rescale=[0.0,0.0,0.0,1.0,1.0,1.0], buff=3):
        """
        transform the dmd modes into vtk multiblock vtk \n
        vars: list of variables to be saved
        block_list: list of block names
        snap_type: type of snapshot
        grid3d: 3D grid
        rescale: rescale the grid points. [x_shift, y_shift, z_shift, x_norm, y_norm, z_norm]
        """

# ----- setup vtk blocks

        vtk_blocks = list()
        i_start    = 0
        
        for bl_num in block_list:
            
            g = grid3d.g[bl_num-1]
            
            px = (g.px[buff:-buff]+rescale[0])/rescale[3]
            py = (g.py[buff:-buff]+rescale[1])/rescale[4]
            pz = (g.pz[buff:-buff]+rescale[2])/rescale[5]
            
            # modify the grid points arrays based on snapshot type
            if   snap_type == 'block': pass
            elif snap_type == 'X':     px = np.array([0.0])
            elif snap_type == 'Z':     pz = np.array([0.0])
            elif snap_type == 'Y' or snap_type == 'W': py = np.array([0.0])
            
            # build one vtk block
            bl_vtk = create_3d_vtkRectilinearGrid( px, py, pz )            
            n_cells = bl_vtk.GetNumberOfCells()
            
# --------- reconstructed data

            step = self.step
            
            for i in range(step):
                
                header = f"recons_{i:05d}"
                
                for j, var in enumerate(vars):
                
                    data_header = header + '_' + var
                    i_s = i_start + j*n_cells
                    i_e = i_start + (j+1)*n_cells
                    #print(i_s, i_e)
                    
                    databuff = self.recons_data.real[i_s:i_e,i]
                    
                    bl_vtk = add_var_vtkRectilinearGrid( bl_vtk, data_header, databuff)

# --------- modes

            # n_modes = len(self.indxes)
            
            # for i in range(n_modes):
                    
            #         header = f"phi_{self.indxes[i]:05d}"
                    
            #         for j, var in enumerate(vars):
                    
            #             data_header = header + '_' + var
            #             i_s = i_start + j*n_cells
            #             i_e = i_start + (j+1)*n_cells
            #             databuff = self.Phis[i_s:i_e,i]
                        
            #             bl_vtk = add_var_vtkRectilinearGrid( bl_vtk, data_header, databuff)

            vtk_blocks.append( bl_vtk )
            i_start += len(vars)*n_cells
        # build the multiple block dataset
        
        dataset = create_multiblock_dataset( vtk_blocks )
        
        return dataset
        


# ----------------------------------------------------------------------
# >>> Match mesh                                                (Nr.)
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
# ----------------------------------------------------------------------

    def match_mesh( self, df_mesh, snap_type ):
        
        header_re = [ f"recons_{i:05d}" for i in range( self.step ) ]
        header_phi = [ f"phi_{indx:05d}" for indx in self.indxes ]
        
        df_recons = pd.DataFrame( self.recons_data, columns=header_re )
        df_phis = pd.DataFrame( self.Phis, columns=header_phi )
        
        # stack two dataframe based on columns
        
        self.df_modes = pd.concat([df_mesh, df_recons, df_phis], axis=1)
        
        self.df_modes['recons_std_dmd'] = self.recons_std_dmd
        
        del self.recons_std_dmd
        del self.recons_data
        del self.Phis
        
        if snap_type == 'block':
            
            self.df_modes.sort_values(by=['z','y','x'],inplace=True)
            self.GX_header = ['x','y','z']
        
        elif snap_type == 'X':
            
            self.df_modes.sort_values(by=['z','y'],inplace=True)
            self.GX_header = ['y','z']

        elif snap_type == 'Y' or snap_type == 'W':
            
            self.df_modes.sort_values(by=['z','x'],inplace=True)
            self.GX_header = ['x','z']
        
        elif snap_type == 'Z':
            
            self.df_modes.sort_values(by=['y','x'],inplace=True)
            self.GX_header = ['x','y']
        
        
        
# ----------------------------------------------------------------------
# >>> Interpolate reconstructed data                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def interp_recons( self, header, snap_type, coord_shift=True ):
        
        # shift the coordinates
        
        x_imp   = float(self.case_parameters.get('x_imp'))
        H       = float(self.case_parameters.get('H'))
        H_md    = float(self.case_parameters.get('H_md'))
        delta_0 = float(self.case_parameters.get('delta_0'))
        
        if coord_shift:
            
            if snap_type == 'X':
                raise ValueError("X slice not supported.")
            
            elif snap_type == 'Y':
                xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
                yy =   np.array(self.grids_interp[1]) / delta_0
            
            elif snap_type == 'Z':
                xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
                yy = ( np.array(self.grids_interp[1]) + H - H_md ) / delta_0                    

            
        else:
            
            if snap_type == 'X':
                raise ValueError("X slice not supported.")
            
            elif snap_type == 'Y':
                xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
                yy =   np.array(self.grids_interp[1]) / delta_0
            
            elif snap_type == 'Z':
                xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
                yy =   np.array(self.grids_interp[1]) / delta_0   
        
        # reconstruct modal variable; interpolate
        
        v = np.array( self.df_modes[ header ] )
        
        # when doing interpolation, do not shift the coordinate avoiding 
        # necessary steps.
        
        x_origin = self.df_modes[self.GX_header[0]]
        y_origin = self.df_modes[self.GX_header[1]]
         
        v =   griddata( (x_origin, y_origin),
                        v.real,
                        (self.grids_interp[0],self.grids_interp[1]),
                        method='linear' )
        
        return xx, yy, v


# ----------------------------------------------------------------------
# >>> Show mode                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/09  - created
#
# Desc
#   1. shift coordinate
#   2. advance in phase
#   3. interpolate onto grids_interp
# ----------------------------------------------------------------------

    def interp_mode( self, indx, snap_type, phase=1.0+0.0j, coord_shift=True ):
        
        # shift the coordinates
        
        x_imp   = float(self.case_parameters.get('x_imp'))
        H       = float(self.case_parameters.get('H'))
        H_md    = float(self.case_parameters.get('H_md'))
        delta_0 = float(self.case_parameters.get('delta_0'))
        
        if coord_shift:
            
            if snap_type == 'X':
                raise ValueError("X slice not supported.")
            
            elif snap_type == 'Y':
                xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
                yy =  np.array(self.grids_interp[1]) / delta_0
            
            elif snap_type == 'Z':
                xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
                yy = ( np.array(self.grids_interp[1]) + H - H_md ) / delta_0                    


        else:
            
            if snap_type == 'X':
                raise ValueError("X slice not supported.")
            
            elif snap_type == 'Y':
                xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
                yy =  np.array(self.grids_interp[1]) / delta_0
            
            elif snap_type == 'Z':
                xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
                yy =   np.array(self.grids_interp[1]) / delta_0  
        
        
        # reconstruct modal variable; interpolate

        header = f"phi_{indx:05d}"
        
        phi = np.array( self.df_modes[ header ] )
        
        alpha_pol = self.df_ind.loc[self.df_ind['indxes']==indx,
                                    'alphas_pol'].iloc[0]
        
        v = phi * alpha_pol * phase
        
        x_origin = self.df_modes[self.GX_header[0]]
        y_origin = self.df_modes[self.GX_header[1]]
        
        v = griddata( (x_origin, y_origin),
                      v.real,
                      (self.grids_interp[0],self.grids_interp[1]),
                      method='linear' )
        
        return xx, yy, v



# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
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
# ----------------------------------------------------------------------

def Testing():

    pass



# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
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
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
