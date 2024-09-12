# -*- coding: utf-8 -*-
'''
@File    :   solution.py
@Time    :   2024/04/30 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import numpy             as     np
import pandas            as     pd
from   typing            import List, Dict, BinaryIO, Any

from   .io_binary        import read_int_bin
from   .io_binary        import read_flt_bin
from   .io_binary        import read_log_bin
from   .io_binary        import write_flt_bin
from   .io_binary        import write_int_bin
from   .io_binary        import write_log_bin
from   .grid             import GridData
from   .block            import SolutionBlock
from   .timer            import timer

class Solution:
    
    def __init__(self, file_path=None, levelset=False):
        
        if file_path is not None:
        
            self._file_path = file_path
            self._fsize     = os.stat( file_path ).st_size
        
        else: 
        
            self._file_path = None
            self._fsize     = None
        
        self._levelset = levelset
        
        # blocks 
        
        self.bl = []
        
        # list of position pointer to var_data chunk start
        
        self.pos_var_start = []
        
        # if verbose
        
        self.verbose = False

# ----------------------------------------------------------------------
# >>> read solution header                                        (Nr.)
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

    def read_header(self, fs: BinaryIO):
        
        sin, slg, sfl, pos = 4,4,8,0
        
# ----- read integers
        
        self._format = read_int_bin(fs.read(sin), sin)
        
        if self._format > 4: raise ValueError('Invalid solution header format!')
        
        # format: 1,2,3,4 -> n_int: 5,2,7,7
        n_int = 7 if self._format > 2 else (5 if self._format > 1 else 2)
        buf_int = read_int_bin(fs.read(n_int*sin), sin)
        
        self._npv    = buf_int[0]
        self._itstep = buf_int[1]
        
        if self._format > 1:
            self._wdamp_included, self._wdamp = True, buf_int[2:5]
        else:
            self._wdamp_included, self._wdamp = False, np.zeros((3), np.int32)
        
        if self._format > 2:
            self._multispecies, self._nscalars, self._nspecies = True, buf_int[5], buf_int[6]
        else:
            self._multispecies, self._nscalars, self._nspecies = False, 0, 1
        
        pos += (n_int + 1)*sin
        
# ----- read floats

        self._time = read_flt_bin(fs.read(sfl), sfl)
        pos += sfl
        
# ----- read logicals

        n_log = 17
        if self._format > 3:
            buf_log = read_log_bin(fs.read(n_log*slg), slg)
            pos += n_log*slg
        else:
            buf_log = np.zeros((n_log), np.bool_)
        self._staggered         = buf_log[ 0]
        self._compressible      = buf_log[ 1]
        self._Boussinesq        = buf_log[ 2]
        self._potential_flow    = buf_log[ 3]
        self._viscous_flow      = buf_log[ 4]
        self._barotropic        = buf_log[ 5]
        self._pressure_equation = buf_log[ 6]
        self._conserve_species  = buf_log[ 7]
        self._conserve_scalars  = buf_log[ 8]
        self._const_thermo      = buf_log[ 9]
        self._real_thermo       = buf_log[10]
        self._breacting         = buf_log[11]
        self._real_chemistry    = buf_log[12]
        self._sensible_enthalpy = buf_log[13]
        self._vle               = buf_log[14]
        self._mutationpp_lib    = buf_log[15]
        self._mpp_state_TTV     = buf_log[16]

        # header size
        
        self._header_size = pos
        
# ----- get the number of variables
        
        if self._compressible: self._var = ['rhou','rhov', 'rhow', 'rho', 'rhoE']
        else:
            self._var = []
            for i in range( self._npv ): self._var.extend([ f"var_{'00'[:-len(str(i+1))]}{i+1}" ])
        
        if self._compressible and self._const_thermo:
            self._n_fields = self._npv + 1  # include p
            self._var.extend(['p'])
            if self._nspecies > 1: 
                self._var.extend(f'Kval_{i+1}' for i in range(self._nspecies)) # Inlcude Kvals
        else: raise ValueError("Not yet supported for processing")
        
        if self._levelset:
            self._n_fields += 1
            self._var.extend( ['levelset'] )
        
        if self._wdamp_included:
            self._n_fields += 1
            self._var.extend( ['wdamp'] )


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
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

    def read_data( self, fs: BinaryIO, grid: GridData, blocklist, sel_vars ):
        
        """
        fs     : opened file object
        grid   : grid data object
        blocklist : list of block numbers
        sel_vars  : list of selected variables
        """
        
        end_of_file = False
        
        # new position after reading header
        
        fs.seek( self._header_size )
        pos = self._header_size
        
        # number of variables in this solution file
        
        n_vars = self._n_fields
        
        # read in block one by one
        # only blocks in fill list will be filled with data chunk
        while not end_of_file:
            
            self.bl.append(SolutionBlock(fs, blocklist, n_vars, sel_vars, grid, self.verbose))
            pos = pos + self.bl[-1].size
            
            if pos >= self._fsize: end_of_file = True
        
        
# ----------------------------------------------------------------------
# >>> print solution information                                 (Nr.)
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

    def _print_info(self):
        
        print(f'\n Solution file info')
        print(f' Path    := ...{self._file_path[min(-40,-len(self._file_path)):]}')
        print(f' Size    := {self._fsize} bytes')
        print(f' Format  := {self._format}')
        print(f' Time    := {self._time:.3f}')
        print(f' Step    := {self._itstep}')
        print(f' Vars Tr := {self._npv}')
        print(f" Scalars := {self._nscalars}")
        print(f" Species := {self._nspecies}")
        print(f" Fields  := {self._var}, {self._n_fields} in total\n")
        
        if self.verbose:
            print(" Solution file flags")
            print(f" - Staggered ...{'yes' if self._staggered else 'no'}")
            print(f" - Compressible ...{'yes' if self._compressible else 'no'}")
            print(f" - Boussinesq ...{'yes' if self._Boussinesq else 'no'}")
            print(f" - Potential flow ...{'yes' if self._potential_flow else 'no'}")
            print(f" - Viscous flow ...{'yes' if self._viscous_flow else 'no'}")
            print(f" - Barotropic ...{'yes' if self._barotropic else 'no'}")
            print(f" - Pressure equation ...{'yes' if self._pressure_equation else 'no'}")
            print(f" - Conserve species ...{'yes' if self._conserve_species else 'no'}")
            print(f" - Conserve scalars ...{'yes' if self._conserve_scalars else 'no'}")
            print(f" - Const thermo ...{'yes' if self._const_thermo else 'no'}")
            print(f" - Real thermo ...{'yes' if self._real_thermo else 'no'}")
            print(f" - Breacting ...{'yes' if self._breacting else 'no'}")
            print(f" - Real chemistry ...{'yes' if self._real_chemistry else 'no'}")
            print(f" - Sensible enthalpy ...{'yes' if self._sensible_enthalpy else 'no'}")
            print(f" - VLE ...{'yes' if self._vle else 'no'}")
            print(f" - Mutationpp lib ...{'yes' if self._mutationpp_lib else 'no'}")
            print(f" - Mpp state TTV ...{'yes' if self._mpp_state_TTV else 'no'}")
            print(f" - Levelset ...{'yes' if self._levelset else 'no'}\n")    


# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/04/30  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    filename = '/media/wencan/Expansion/temp/240210/results/results.bin'
    gfile = '/media/wencan/Expansion/temp/240210/results/inca_grid.bin'
    
    bbox = [-200,200,-2,100,-12,12]
    
    grid = GridData(gfile)
    grid.read_grid()
    blocklist = grid.select_blockgrids(bbox)
    
    sol = Solution(filename)
    
    with open(filename, 'rb') as fs:
        sol.read_header(fs)
        sol._print_info()
        sol.read_data(fs,grid,blocklist,['rhou','rhow'])




# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/04/30  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
