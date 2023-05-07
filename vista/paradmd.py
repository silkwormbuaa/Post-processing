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
        
        self.snap_struct_file = snap_dir + '/snap_struct.csv'
        
        self.snap_info_file = snap_dir + '/snap_info.dat'
        
        
        if self.rank == 0:
            
            if not (os.path.exists( self.snap_info_file ) and 
                    os.path.exists( self.snap_struct_file )):
                
                raise FileNotFoundError("no paradmd pre-processing files")
        
            else: print("Checked pre-processing files.")
            
        
        # some parameters need to be defined first
        
        # - snap_info variables
        
        self.kind = None
        self.n_bl = None
        self.snap_type = None
        self.slic_type = None
        self.vars_name = None
        
        # - snap_struct lists
        
        self.bl_num = None
        self.pos = None
        self.size = None
        self.N1 = None
        self.N2 = None
        self.N3 = None
        

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
        
        # Root read the info_file
        
        if self.rank == 0:
        
            with open( self.snap_info_file ) as fi:
                
                self.kind = int(fi.readline().strip().split()[1])
                
                self.n_bl = int(fi.readline().strip().split()[1])
                
                self.snap_type = fi.readline().strip().split()[1]
                
                if self.snap_type == 'slice':
                    self.slic_type = fi.readline().strip().split()[1]
                
                self.vars_name = fi.readline().strip().split()[1:]
            
        
        # Root broadcast the info to other processes
        
        self.kind = self.comm.bcast( self.kind, root=0 )
        
        self.n_bl = self.comm.bcast( self.n_bl, root=0 )
        
        self.snap_type = self.comm.bcast( self.snap_type, root=0 )
        
        self.slic_type = self.comm.bcast( self.slic_type, root=0 )
        
        self.vars_name = self.comm.bcast( self.vars_name, root=0 )
                
                
            
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
        
        # Read snap_struct.dat file into a pandas dataframe
        
        if self.rank == 0:
            
            df = pd.read_csv( self.snap_struct_file, sep=' ' )
            
            self.bl_num = np.array( df['bl_num'] )
            
            self.pos = np.array( df['pos'] )
            
            self.size = np.array( df['size'] )
            
            self.N1 = np.array( df['N1'] )
            
            self.N2 = np.array( df['N2'] )
            
            self.N3 = np.array( df['N3'] )
        
        
        # Broadcast to other processors. 
        # (Can also use single-segment buffer interface with Bcast)
            
        self.bl_num = self.comm.bcast( self.bl_num, root=0 )
        
        self.pos    = self.comm.bcast( self.pos, root=0 )
        
        self.size   = self.comm.bcast( self.size, root=0 )
        
        self.N1     = self.comm.bcast( self.N1, root=0 )
        
        self.N2     = self.comm.bcast( self.N2, root=0 )
        
        self.N3     = self.comm.bcast( self.N3, root=0 )
    

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
    
    paradmd.comm.barrier()
    
    paradmd.read_info_file()
    
    paradmd.read_struct_file()
    
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.bl_num[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.pos[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.size[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.N1[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.N2[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.N3[-1]}')
    
'''
    print(f'This is proc {paradmd.rank}, I got the kind is {paradmd.kind}.')
    print(f'This is proc {paradmd.rank}, I got the kind is {paradmd.n_bl}.')
    print(f'This is proc {paradmd.rank}, I got the kind is {paradmd.snap_type}.')
    print(f'This is proc {paradmd.rank}, I got the kind is {paradmd.slic_type}.')
    print(f'This is proc {paradmd.rank}, I got the kind is {paradmd.vars_name}.')
'''    
                                                     



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

    with timer('Testing '):
        
        Testing()
