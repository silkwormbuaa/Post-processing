#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   params.py
@Time    :   2024/10/23 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   define params class and related functions
'''


import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )


class Params:

# ----------------------------------------------------------------------
# >>> Initialize Params class                                    (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def __init__( self, paramsfile ):
        
        """
        Initiate the class with the case_parameters file
        
        return: self.params, dictionary of parameters
        """

        self.paramsfile = paramsfile
        
        self.params = self.read_params()
        

# ----------------------------------------------------------------------
# >>> read params                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_params( self ):
        
        # create a dictionary for parameter pairs
        
        parameters = {}

        if not os.path.exists( self.paramsfile ):
            raise FileExistsError("Please set case_parameters file!")

        # read file content line by line
        
        with open( self.paramsfile, 'r') as f:
            
            lines = f.readlines()

            for line in lines:
                
                # ignore notation line and space line
                
                if line.strip() and not line.startswith("#"):
                    
                    # split from the '#'
                    
                    line = line.split("#")[0]
                    
                    # split from the first '='
                    
                    key, value = line.strip().split("=",1)
                    
                    # remove extra space and comma
                    
                    key = key.strip()
                    value = value.strip().rstrip(",")
                    
                    # add parameter pair into the dictionary
                    
                    parameters[key] = value
        
        return parameters

# ----------------------------------------------------------------------
# >>> define each attribute                                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

# ------------------- Basic case info ---------------------------------

    @property
    def casecode( self ):
        return self.params.get('casecode')
    
    @property
    def tag( self ):
        return self.params.get('tag')
    
    @property
    def roughwall( self ):
        return True if str(self.params.get('roughwall')).lower() == 'true' else False


# -------------------- Flow setup parameters --------------------------

    @property
    def u_ref( self ):
        return float( self.params.get('u_ref') )
    
    @property
    def p_ref( self ):
        return float( self.params.get('p_ref') )
    
    @property
    def rho_ref( self ):
        return float( self.params.get('rho_ref') )
    
    @property
    def T_ref( self ):
        return float( self.params.get('T_ref') )
    
    @property
    def Re_ref( self ):
        return float( self.params.get('Re_ref') )
    
    @property
    def delta_0( self ):
        return float( self.params.get('delta_0') )
    
    @property
    def visc_law( self ):
        return self.params.get('visc_law')


# ------------------- probe setup parameters --------------------------

    @property
    def prb_ridge_index( self ):
        i_s, i_e = self.params.get('prb_ridge_index').split(',')
        return int(i_s), int(i_e)

    @property
    def prb_valley_index( self ):
        i_s, i_e = self.params.get('prb_valley_index').split(',')
        return int(i_s), int(i_e)

    @property
    def prb_upspan_index( self ):
        i_s, i_e = self.params.get('prb_upspan_index').split(',')
        return int(i_s), int(i_e)
        
    @property
    def prb_downspan_index( self ):
        i_s, i_e = self.params.get('prb_downspan_index').split(',')
        return int(i_s), int(i_e)
    
    @property
    def prb_upvline_index( self ):
        i_s, i_e = self.params.get('prb_upvline_index').split(',')
        return int(i_s), int(i_e)
    
    @property
    def prb_downvline_index( self ):
        i_s, i_e = self.params.get('prb_downvline_index').split(',')
        return int(i_s), int(i_e)

        
# ------------------- Numerical setup parameters ----------------------

    @property
    def dt_snap( self ):
        return float( self.params.get('dt_snap') )


# ---------------------- Geometry parameters --------------------------
    
    @property
    def H( self ):
        return float( self.params.get('H') )
    
    @property
    def H_md( self ):
        return float( self.params.get('H_md') )
    
    @property
    def x_imp( self ):
        return float( self.params.get('x_imp') )
    
    @property
    def period( self ):
        return float( self.params.get('period') )
    
    @property
    def D( self ):
        return float( self.params.get('D') )
    

# --------------- case specific flow conditions -----------------------
    
    @property
    def x_sep( self ):
        """in shifted and normalized coordinate"""
        return float( self.params.get('x_sep') )
    
    @property
    def x_att( self ):
        """in shifted and normalized coordinate"""
        return float( self.params.get('x_att') )
    
    @property
    def lsep( self ):
        """in shifted and normalized coordinate"""
        return self.x_att - self.x_sep
    
    @property
    def Lsep( self ):
        """in physical unit"""
        return self.lsep * self.delta_0

    @property
    def x_pfmax( self ):
        """in shifted and normalized coordinate"""
        return float( self.params.get('x_pfmax') )
    

# ------------------- INCA version led difference ------------------------
    
    @property
    def prb_withT( self ):
        return True if self.params.get('prb_withT').lower() == 'true' else False


# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    paramfile = '/media/wencanwu/Seagate Expansion Drive1/temp/231124/supplements/case_parameters'

    params = Params( paramfile )
    
    print( params.prb_ridge_index )


# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
