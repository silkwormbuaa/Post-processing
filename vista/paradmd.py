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

import gc 

import numpy             as np

import pandas            as pd

from   mpi4py            import MPI

sys.path.append('..')

from   utils.timer       import timer

from   utils.tools       import get_filelist

from   utils.init_empty  import init_1Dflt_empty

from   utils.init_empty  import init_2Dflt_empty

from   utils.init_empty  import init_2Dcmx_empty

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
        self.type = None
        self.slic_type = None
        self.vars_name = None
        self.n_var = None
        
        # - snap_struct lists
        
        self.bl_num = None
        self.pos = None
        self.size = None
        self.N1 = None
        self.N2 = None
        self.N3 = None
        
        # - select which variable
        
        self.select = None
        
        # - snapshots matrix for dmd
        
        self.snapshots = []
        
        self.Ns = None
        


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
                
                self.type = fi.readline().strip().split()[1]
                
                if self.type == 'slice':
                    self.slic_type = fi.readline().strip().split()[1]
                
                self.vars_name = fi.readline().strip().split()[1:]
                
                self.n_var = len( self.vars_name )
            
        
        # Root broadcast the info to other processes
        
        self.kind = self.comm.bcast( self.kind, root=0 )
        
        self.n_bl = self.comm.bcast( self.n_bl, root=0 )
        
        self.type = self.comm.bcast( self.type, root=0 )
        
        self.slic_type = self.comm.bcast( self.slic_type, root=0 )
        
        self.vars_name = self.comm.bcast( self.vars_name, root=0 )
                
        self.n_var     = self.comm.bcast( self.n_var, root=0 )
        
        
            
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
# - Calculate number of blocks a proc should take care: n_bl_local 
#
# - Calculate start index 'i_start' within 0,...,n_bl-1
#   or i_start = None and n_bl_local = 0 if there is no block
#   (never should happen)
#
# ----------------------------------------------------------------------

    def assign_block( self ):
        
        # algorithm divides number of blocks in a dealing way
        
        n_bl_local = self.n_bl // self.n_procs # floor dividing
        
        left_procs = self.n_bl - n_bl_local*self.n_procs
        

        if ( self.rank < left_procs ):
            
            n_bl_local = n_bl_local + 1
            
            i_start = 0 + self.rank*n_bl_local
            
            i_end = i_start + n_bl_local

            
        elif ( n_bl_local > 0 ):
            
            i_start = 0 + self.rank*n_bl_local + left_procs
            
            i_end = i_start + n_bl_local


        else:
            
            n_bl_local = 0
            
            i_start = None
            
            i_end = None
        
        
        # Passing out the parameters
        
        self.n_bl_local = n_bl_local
        
        self.i_start = i_start
        
        self.i_end = i_end
        
        # For the last proc's i_end, i_end = n_bl is out of bound 
        # if it is used to indexing or slicing out a block.
        # But luckily we won't use it. :)
            
    

# ----------------------------------------------------------------------
# >>> parallel read data chunk                                    (Nr.)
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
#   - form and return vector
#
# ----------------------------------------------------------------------

    def para_read_data( self, filename ):
        
        # MPI-IO parallel reading different part of a snapshot
        
        fh = MPI.File.Open( self.comm, filename, amode=MPI.MODE_RDONLY )
        
        data = None
        
        # Set views for each blocks
        # etype and filetype are byte!(can customize, but not necessary

        fh.Set_view( self.pos[self.i_start], MPI.BYTE, MPI.BYTE )
        
        
        # Read data from each block
        
        for bl in range( self.i_start, self.i_end ):
            
            # Move pointer first
            
            pos_move = self.pos[bl] - self.pos[self.i_start]
            
            fh.Seek( pos_move, whence = MPI.SEEK_SET )
            
            
            # Create a buffer and read data chunk into buffer
            buff_data = init_1Dflt_empty( int(self.size[bl]/self.kind),
                                          self.kind )
#            buff_data = np.empty( int(self.size[bl]/self.kind), dtype=np.float32)
        
            fh.Read( buff_data )
            
            
            # Drop ghost cells, other variables, and stack along blocks
            
            if data is None: 
                
                data = self.drop_ghost( buff_data, bl )
#                print(f'data shape is {np.shape(data)}')
                
            else: 
                
                data = np.hstack((data, self.drop_ghost(buff_data, bl)))
        
            del buff_data
#            gc.collect()
        
        self.snapshots.append( data )

        fh.Close()
        

# ----------------------------------------------------------------------
# >>> Drop ghost                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/07  - created
#
# Desc
#
# 
# ----------------------------------------------------------------------

    def drop_ghost( self, buff_data, bl, ghost=3 ):
        
        # Get the shape of current data chunk
        
        N1 = self.N1[bl]
        
        N2 = self.N2[bl]
        
        N3 = self.N3[bl]
        
        
        if self.type == 'block':
            
            Nx = N1 - ghost*3
            Ny = N2 - ghost*3
            Nz = N3 - ghost*3
        
            buff_data = buff_data.reshape(( self.n_var, N1, N2, N3 ))
            
            buff_data = buff_data[ : , ghost:-ghost, ghost:-ghost, ghost:-ghost]
        
        
        elif self.type == 'slice':
            
            if self.slic_type == 'X':
                
                Nx = 1
                Ny = N2 - ghost*2
                Nz = N3 - ghost*2
                
                buff_data = buff_data.reshape(( self.n_var, N2, N3 ))
                
                buff_data = buff_data[ : , ghost:-ghost, ghost:-ghost]
                

            if self.slic_type == 'Y' or self.slic_type == 'W':
                
                Nx = N1 - ghost*2
                Ny = 1
                Nz = N3 - ghost*2
                
                buff_data = buff_data.reshape(( self.n_var, N1, N3 ))
                
                buff_data = buff_data[ : , ghost:-ghost, ghost:-ghost]
                

            if self.slic_type == 'Z':
                
                Nx = N1 - ghost*2
                Ny = N2 - ghost*2
                Nz = 1
                
                buff_data = buff_data.reshape(( self.n_var, N1, N2 ))
                
                buff_data = buff_data[ : , ghost:-ghost, ghost:-ghost]
        
        
        # reshape to long data vectors of different variable
        
        buff_data = buff_data.reshape(( self.n_var, Nx*Ny*Nz ))

        if self.select == 'u':    buff_data = buff_data[0,:]
        elif self.select == 'v':  buff_data = buff_data[1,:]
        elif self.select == 'w':  buff_data = buff_data[2,:]
        elif self.select == 'T':  buff_data = buff_data[3,:]
        elif self.select == 'p':  buff_data = buff_data[4,:]
        elif self.select == 'cf': buff_data = buff_data[5,:]
        else: raise ValueError('The selected variable does not exist.')
        
        return buff_data
        


# ----------------------------------------------------------------------
# >>> Parallel DMD                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def do_paradmd( self ):
        
        # Check if the snapshots data are available
        
        self.N_t = len( self.snapshots )
        
        if self.Ns == 0:
            
            raise ValueError('Please read snapshots first!')     
        
        
        # just for shortening name
        
        Ns = self.N_t - 1
        n_procs = self.n_procs
        rank = self.rank
        
        # ==============================================================
        # Do QR factorization in parallel
        # ==============================================================
        
        # 1. Direct QR factorization of Ai(1 => Ns)
                
        Ai = np.array(self.snapshots).T
        
        del self.snapshots
        
        Q1i, Ri = np.linalg.qr( Ai[:,:-1], mode='reduced' )
        
        
        # 2. Build the buffer also form the Rprime matrix
        
        Rprime =  init_1Dflt_empty( Ns*Ns*n_procs, self.kind )

        Rprime[ Ns*Ns*rank : Ns*Ns*(rank+1) ] = Ri.ravel()
        
        for r in range( n_procs ):
            
            self.comm.Bcast( Rprime[ Ns*Ns*r : Ns*Ns*(r+1) ], root=r )
        
        del Ri
        gc.collect()
        
        
        # 3. QR factorization of Rprime and get R, Q2i
        
        Q2i, R = np.linalg.qr( Rprime.reshape(Ns*n_procs, Ns), mode='reduced')

        del Rprime
        gc.collect()
        
        
        # 4. Compute Qi = Q1i x Q2i
        
        Qi = np.matmul( Q1i, Q2i[Ns*rank:Ns*(rank+1), : ] )
        
        del Q1i, Q2i
        gc.collect()
        
        
        # ==============================================================
        # Compute SVD of Ai = Qi x R
        # ==============================================================
        
        # 1. Compute SVD of R
        
        Ur, S, Vt = np.linalg.svd( R, full_matrices=False )
        
        
        # 2. Build up sigma
        
        sigma = np.diag( S )
        
        
        # 3. Compute Ui = Qi x Ur
        
        Ui = np.matmul( Qi, Ur )
        
        del Ur, S
        gc.collect()
        
        
        # ==============================================================
        # Solve eigenvalue problem of S tilde
        # ==============================================================
        
        # 1. Project S onto POD basis of Ai(2 => Ns+1) 
        #    Ui* x Ai(2 => Ns+1)
        
        B = np.matmul( Ui.conj().T, Ai[:,1:] )
        
        
        # 2. Sum Ui* x Ai(2 => Ns+1) over all processors
        
        C = init_2Dflt_empty( Ns, Ns, self.kind )
                
        self.comm.Allreduce( B, C, op=MPI.SUM )
        
        
        # 3. Compute C x V = Ui* x Ai(2 => Ns+1) x V
        
        '''why C.transpose()'''
        
        B = np.array( np.matmul( C, Vt.conj().T ), order='F' )
        
        del C
        gc.collect()
        
        
        # 4. Compute S tilde
        
        Si_td = np.matmul( B, np.linalg.inv( sigma) )
        
        del B 
        gc.collect()
        
        
        # 5. Compute the eig(Si_st), with Si_d = Ui* x Ai(2=>Ns) x V x sigma^-1
        
        mu, Y = np.linalg.eig( Si_td )
        
        del Ai
        gc.collect()
        
#        print(mu)
        
        # ==============================================================
        # DMD modes in real space and amplitude
        # ==============================================================
        
        # DMD modes
        
        Phi_i = np.matmul( Ui, Y )
        
        
        # Build Vandermonde matrix
        
        Vand = init_2Dcmx_empty( Ns, Ns, self.kind*2 )
        
        Vand[:,0] = 1.0
        
        Vand[:,1] = mu

        for j in range(2, Ns): 
            Vand[:,j] = Vand[:,j-1] * mu
            
        
        # Build objective function matrices
        P = np.multiply( np.dot( Y.conj().T, Y )
                    , np.conj( np.dot( Vand, Vand.conj().T ) ) )
        q = np.conj( np.diag( np.dot( np.dot( np.dot( Vand, Vt.conj().T ), 
                                             sigma.conj().T ), Y ) ) )
        s = np.trace( np.dot( sigma.conj().T, sigma ) )

        
        # Find amplitudes directly
        alphas = np.linalg.solve( P, q )
        
        
        self.mu = mu
        
        self.alphas = np.abs( alphas ) 
        
        self.freq = np.angle( mu )/0.01/(2*np.pi)
        
        
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
    
    snap_file = snap_dir+'/snap_test/snapshot_00600031/snapshot_W_002.bin'
    
    paradmd = ParaDmd( snap_dir )
    
    paradmd.comm.barrier()
    
    paradmd.read_info_file()
    
    paradmd.read_struct_file()
    
    paradmd.assign_block()
    
    print('I am process %d out of %d, got the total block number %d.\
            I take care of blocks from %d(%d) .. %d(%d)'
            %(paradmd.rank, paradmd.n_procs, paradmd.n_bl,
            paradmd.i_start+1, paradmd.bl_num[paradmd.i_start], 
            paradmd.i_end, paradmd.bl_num[paradmd.i_end-1])) 
    
    paradmd.comm.barrier()
    
    snap_files = None
    
    if paradmd.rank == 0:
        
        snap_files = get_filelist( snap_dir + '/snap_test' )
        
    snap_files = paradmd.comm.bcast( snap_files, root=0 )
    
    paradmd.select = 'p'
    
    for snap_file in snap_files:
        
        paradmd.para_read_data( snap_file )
    
    print(np.shape( paradmd.snapshots ))
'''
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.bl_num[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.pos[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.size[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.N1[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.N2[-1]}')
    print(f'This is proc {paradmd.rank}, I got bl_num {paradmd.N3[-1]}')
    
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
