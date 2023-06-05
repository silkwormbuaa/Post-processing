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

    def __init__( self ):
        
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
        self.Phi      = None
        


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

    def add_modes( self, mode ):
        
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

        Phis = [mode.Phi for mode in self.modes]

        alpha_pols = [mode.alpha_pol for mode in self.modes]
        
        mus = [mode.mu for mode in self.modes] 
        
        vand = np.vander( mus, step, increasing=True )
        
        if not self.modes:
        
            raise ValueError("No mode is available!")
        
        else:
            
            self.recons_data = Phis.T @ np.diag(alpha_pols) @ vand



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
