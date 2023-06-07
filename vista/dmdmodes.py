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
        indx:        index of the mode among all standard DMD modes
        alpha:      amplitude from standard DMD
        alpha_sp:   amplitude before polishing from SPDMD 
        alpha_pol:  amplitude after polishing from SPDMD
        mu:          eigenvalues from standard DMD
        '''
        self.indx      = None
        self.alpha     = None
        self.alpha_sp  = None
        self.alpha_pol = None
        self.mu        = None
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

        self.recons_data = None


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
        
        alphas = [mode.alpha for mode in self.modes]

        alpha_pols = [mode.alpha_pol for mode in self.modes]
        
        mus = [mode.mu for mode in self.modes] 
        
        vand = np.vander( mus, self.step, increasing=True )

        
        if not self.modes:
        
            raise ValueError("No mode is available!")
        
        else:
            
            self.recons_data = self.Phis @ np.diag(alphas) @ vand



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

    def match_mesh( self, df, snap_type ):
        
        header_re = [ f"recons_{i:05d}" for i in range( self.step )]
        header_phi = [ f"phi_{indx:05d}" for indx in self.indxes]
        
        df_recons = pd.DataFrame( self.recons_data, columns=header_re )
        df_phis = pd.DataFrame( self.Phis, columns=header_phi )
        
        # stack two dataframe based on columns
        
        self.df_modes = pd.concat([df, df_recons, df_phis], axis=1)
        
        if snap_type == 'block':
            
            self.df_modes.sort_values(by=['z','y','x'],inplace=True)
        
        elif snap_type == 'X':
            
            self.df_modes.sort_values(by=['z','y'],inplace=True)

        elif snap_type == 'Y' or snap_type == 'W':
            
            self.df_modes.sort_values(by=['z','x'],inplace=True)
        
        elif snap_type == 'Z':
            
            self.df_modes.sort_values(by=['y','x'],inplace=True)


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
