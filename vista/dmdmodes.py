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

    def interp_recons( self, header, coord_shift=True ):
        
        # shift the coordinates
        
        if coord_shift:
            
            x_imp   = float(self.case_parameters.get('x_imp'))
            H       = float(self.case_parameters.get('H'))
            H_md    = float(self.case_parameters.get('H_md'))
            delta_0 = float(self.case_parameters.get('delta_0'))
            
            xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
            yy = ( np.array(self.grids_interp[1]) + H - H_md ) / delta_0
            
        else:
            
            xx = np.array( self.grids_interp[0] )
            yy = np.array( self.grids_interp[1] )
        
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

    def interp_mode( self, indx, phase=1.0+0.0j, coord_shift=True ):
        
        # shift the coordinates
        
        if coord_shift:
            
            x_imp   = float(self.case_parameters.get('x_imp'))
            H       = float(self.case_parameters.get('H'))
            H_md    = float(self.case_parameters.get('H_md'))
            delta_0 = float(self.case_parameters.get('delta_0'))
            
            xx = ( np.array(self.grids_interp[0]) - x_imp ) / delta_0
            yy = ( np.array(self.grids_interp[1]) + H - H_md ) / delta_0
            
        else:
            
            xx = np.array( self.grids_interp[0] )
            yy = np.array( self.grids_interp[1] )
        
        
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
