# -*- coding: utf-8 -*-
'''
@File    :   paradmd.py
@Time    :   2023/05/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os

import sys

import numpy             as np

import pandas            as pd

from   mpi4py            import MPI

sys.path.append('..')

from   utils.timer       import timer


class ParaDmd:
# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def __init__( self, snap_dir ):
        
        # Initialize MPI environment
        
        self.comm    = MPI.COMM_WORLD
        self.rank    = self.comm.Get_rank()
        self.n_procs = self.comm.Get_size()
        
        
        # Check if pre-processing files exists
        
        self.snap_dir = snap_dir
        
        snap_struct_file = snap_dir + '/snap_struct.csv'
        
        snap_info_file = snap_dir + '/snap_info.dat'
        
        
        if self.rank == 0:
            
            if not (os.path.exists( snap_info_file ) and 
                    os.path.exists( snap_struct_file )):
                
                raise FileNotFoundError("no paradmd pre-processing files")
        
            else: print("Checked pre-processing files.")
        
        

# ----------------------------------------------------------------------
# >>> Read info file and broadcast                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/06  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_info_file( self ):
        
        pass
    

# ----------------------------------------------------------------------
# >>> Read struct file and broadcast                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/06  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_struct_file( self ):
        
        pass
    

# ----------------------------------------------------------------------
# >>> Assign Works                                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/06  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def assign_block( self ):
        
        pass 
    

# ----------------------------------------------------------------------
# >>> parallel reading block                                    (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/06  - created
#
# Desc
#   - read block data chunk
#   - drop ghost point
#   - form matrix
#
# ----------------------------------------------------------------------



# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    snap_dir = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots'
    
    paradmd = ParaDmd( snap_dir )
                                                     



# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
