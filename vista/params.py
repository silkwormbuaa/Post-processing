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
import re
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
    
    @property
    def p_dyn( self ):
        return 0.5 * self.rho_ref * self.u_ref**2

        
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
        return int( self.params.get('period') )
    
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
    
    @property
    def prb_sep( self ):
        """file index of the probe data at separation point"""
        return int( self.params.get('prb_sep') )
        
    @property
    def prb_att( self ):
        """file index of the probe data at attachment point"""
        return int( self.params.get('prb_att') )
    
    @property
    def prb_pfmax( self ):
        """file index of the probe data at peak pressure point"""
        return int( self.params.get('prb_pfmax') )

# ------------------- INCA version led difference ------------------------
    
    @property
    def prb_withT( self ):
        return True if self.params.get('prb_withT').lower() == 'true' else False


# ------------------- probe setup parameters --------------------------

    def __indexstr_to_list__( indexstr ):
        
        def parse_range(range_str):
            """parse str like '100..105' """
            if ".." in range_str:
                start, end = map(int, range_str.split(".."))
            else:
                return []
            return list(range(start, end + 1))

        def parse_list(list_str):
            """parse str like [20, 30, 40] """
            # remove the square brackets and split by comma
            return list(map(int, list_str.strip("[]").split(", ")))
        
        merged_list = []
        
        # find all the parts with a pair of square brackets
        parts = re.findall(r'\[.*?\]', indexstr)

        for part in parts:
            
            if ".." in part:
                merged_list.extend( parse_range(part.strip("[]")) )
            elif '[]' in part:
                merged_list.extend([])
            else:
                merged_list.extend( parse_list(part.strip("[]")) )
        
        return sorted(merged_list)


    @property
    def prb_ridge_index( self ):
        return self.__indexstr_to_list__( self.params.get('prb_ridge_index') )
    
    @property
    def prb_valley_index( self ):
        return self.__indexstr_to_list__( self.params.get('prb_valley_index') )

    @property
    def prb_upspan_index( self ):
        return self.__indexstr_to_list__( self.params.get('prb_upspan_index') )
        
    @property
    def prb_downspan_index( self ):
        return self.__indexstr_to_list__( self.params.get('prb_downspan_index') )
    
    @property
    def prb_upvline_index( self ):
        return self.__indexstr_to_list__( self.params.get('prb_upvline_index') )
    
    @property
    def prb_downvline_index( self ):
        return self.__indexstr_to_list__( self.params.get('prb_downvline_index') )



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
