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
import gc
import sys
import time
import pickle
import numpy             as     np
import numpy.linalg      as     linalg
import pandas            as     pd
from   mpi4py            import MPI

from   .timer            import timer
from   .colors           import colors   as col
from   .tools            import get_filelist
from   .tools            import to_dictionary
from   .tools            import distribute_mpi_work
from   .init_empty       import init_1Dflt_empty
from   .init_empty       import init_1Dcmx_empty
from   .init_empty       import init_2Dflt_empty
from   .init_empty       import init_2Dcmx_empty
from   .directories      import create_folder

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
        
        """
        snap_dir : directory to '.../snapshots'
        """
        
        # Initialize MPI environment
        
        self.comm    = MPI.COMM_WORLD
        self.rank    = self.comm.Get_rank()
        self.n_procs = self.comm.Get_size()
        
        # Check if pre-processing files exists
        
        self.snap_dir = snap_dir
        
        dmddir = create_folder( snap_dir + '/../postprocess/dmd' )
        
        self.snap_struct_file  = dmddir + '/snap_struct.csv'
        
        self.snap_info_file    = dmddir + '/snap_info.dat'
        
        self.Pqs_file          = dmddir + '/Pqs.pkl'
        
        self.ind_spmode_file   = dmddir + '/ind_spmode.csv'
        
        self.spdmd_result_file = dmddir + '/spdmd_result.pkl'
        
        self.dmdmodes_dir      = dmddir + '/dmdmodes'
        
        
        # some parameters need to be defined first
        
        # - snap_info variables
        
        self.kind = 4
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
        
        self.var_norms = None
        
        # - snapshots matrix for dmd, total number and length of snapshots
        
        self.snapshots = []
        
        self.N_t = None
        
        self.len_snap_local = None
        
        # - time interval
        
        self.dt = None
        
        # - spdmd objective function components
        
        self.P = None
        self.q = None
        self.s = None
        
        # list of output ind_spmode table
        
        self.St = None
        self.beta = None
        self.psi = None
        self.psi_pol = None


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
#   - read info file and broadcast
#   - set default variable values for normalization
#
# ----------------------------------------------------------------------

    def read_info_file( self, info_file=None ):
        
        """
        info_file: If not given, use default self.snap_info_file
        """
        
        # initialize info_file. 
        
        if info_file is None:
            info_file = self.snap_info_file
        
        # Check if files are available.
        
        if self.rank == 0:
            
            if not ( os.path.exists( info_file ) ):
                
                raise FileNotFoundError(col.fg.red,
                                        "no paradmd pre-processing info file",
                                        col.reset)
        
            else: print("Checked pre-processing info file.\n")
        
        
        # Root read the info_file
        
        if self.rank == 0:
        
            with open( info_file ) as fi:
                
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
        
        
        # Set default values 1.0 of normalization
        
        self.var_norms = to_dictionary( self.vars_name, np.ones(self.n_var) )
        
            
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

    def read_struct_file( self, struct_file=None ):
        
        """
        If not given, use default self.snap_info_file
        """
        
        # initialize info_file. 
        
        if struct_file is None:
            struct_file = self.snap_struct_file
        
        # Check if files are available.
        
        if self.rank == 0:
            
            if not ( os.path.exists( struct_file ) ):
                
                raise FileNotFoundError(col.fg.red,
                                        "no paradmd pre-processing struct file",
                                        col.reset)
        
            else: print("Checked pre-processing struct file.\n")        


        # Read snap_struct.dat file into a pandas dataframe
        
        if self.rank == 0:
            
            df = pd.read_csv( struct_file, sep=' ' )
            
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
        
        """
        No parameter needed.
        """
        
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
        
        """
        filename: path to snapshot file
        """
        
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
                
                # initialize data
                
                data = self.drop_ghost( buff_data, bl )
                                
            else: 
                
                data = np.hstack((data, self.drop_ghost(buff_data, bl)))
        
            del buff_data
#            gc.collect()
        
        self.snapshots.append( data )

        fh.Close()
        

# ----------------------------------------------------------------------
# >>> Drop ghost cells and normalize data                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/07  - created
#
# Desc
#   - drop ghost cells of one block
#   - normalize the selected variables
# ----------------------------------------------------------------------

    def drop_ghost( self, buff_data, bl, ghost=3 ):
        
        """
        buff_data: data chunk with ghost cells (n_var*N3*N2*N1 long vector)
        bl: index of current block
        """
        
        # Get the shape of current data chunk
        
        N1 = self.N1[bl]
        
        N2 = self.N2[bl]
        
        N3 = self.N3[bl]
        
        
        if self.type == 'block':
            
            Nx = N1 - ghost*3
            Ny = N2 - ghost*3
            Nz = N3 - ghost*3
        
            buff_data = buff_data.reshape(( self.n_var, N3, N2, N1 ))
            
            buff_data = buff_data[ : , ghost:-ghost, ghost:-ghost, ghost:-ghost]
        
        
        elif self.type == 'slice':
            
            if self.slic_type == 'X':
                
                Nx = 1
                Ny = N2 - ghost*2
                Nz = N3 - ghost*2
                
                buff_data = buff_data.reshape(( self.n_var, N3, N2 ))
                
                buff_data = buff_data[ : , ghost:-ghost, ghost:-ghost]
                

            if self.slic_type == 'Y' or self.slic_type == 'W':
                
                Nx = N1 - ghost*2
                Ny = 1
                Nz = N3 - ghost*2
                
                buff_data = buff_data.reshape(( self.n_var, N3, N1 ))
                
                buff_data = buff_data[ : , ghost:-ghost, ghost:-ghost]
                

            if self.slic_type == 'Z':
                
                Nx = N1 - ghost*2
                Ny = N2 - ghost*2
                Nz = 1
                
                buff_data = buff_data.reshape(( self.n_var, N2, N1 ))
                
                buff_data = buff_data[ : , ghost:-ghost, ghost:-ghost]
        
        
        # Reshape to long data vectors of different variable and 
        # normalize the variables
        
        buff_data = buff_data.reshape(( self.n_var, Nx*Ny*Nz ))

        buff_data = self.select_var( buff_data )
        
        return buff_data
        

# ----------------------------------------------------------------------
# >>> select variables                                           (Nr.)
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

    def select_var( self, buff_data ):

        """
        return: a 1D list of selected variables
        """
        
        temp = []
        vars = self.select
        
        if 'u' in vars:    
            temp.append( buff_data[0,:]/self.var_norms.get('u') )
            
        if 'v' in vars:  
            temp.append( buff_data[1,:]/self.var_norms.get('v') )
            
        if 'w' in vars:  
            temp.append( buff_data[2,:]/self.var_norms.get('w') )
            
        if 'T' in vars:  
            temp.append( buff_data[3,:]/self.var_norms.get('T') )
            
        if 'p' in vars:  
            temp.append( buff_data[4,:]/self.var_norms.get('p') )
            
        if 'cf' in vars: 
            temp.append( buff_data[5,:]/self.var_norms.get('cf') )

        
        return np.array( temp ).ravel()
        

# ----------------------------------------------------------------------
# >>> Non-blocking data exchange                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/06/11  - created
#
# Desc
#    - mpi non-blocking snapshot data exchange 
#    - data are read from each snapshot which composes a long vector,
#      each processor may read different number of snapshots, but in a card
#      shuffling way, so that each processor will have nearly equal number of 
#      snapshots.
# 
#    - after reading, the snapshots are stacked together like below:
#      (n,m) n-number of snapshots, m-length of each snapshot data vector  
#   
#              CPU1           CPU2         CPU3
#
#           s1  s2  s3      s4  s5       s6  s7
#            |   |   |       |   |        |   |
#            |   |   |       |   |        |   |
#            |   |   |       |   |        |   |
#
#    - after the data exchange, each processor will have part of all the snapshots. 
#
#            s1  s2  s3  s4  s5  s6  s7
#    
#      CPU1   |   |   |   |   |   |   |
#    
#      CPU2   |   |   |   |   |   |   |
#    
#      CPU3   |   |   |   |   |   |   |
#
# ----------------------------------------------------------------------

    def non_blocking_data_exchange( self, data ):
        
        comm = self.comm
        rank = self.rank
        size = self.n_procs

        def debug_print(message, flush=False):
            if flush:
                print(f"[Rank {rank}] {message}")
                sys.stdout.flush()
            else:pass
        
        send_requests = []
        recv_requests = []
        
        # check if the input data is 2d np.ndarray
        if not isinstance(data, np.ndarray) or data.ndim != 2:
            raise ValueError("Input data before data exchange must be a 2D numpy array.")

        # local size of the data
        n_cols, n_rows = data.shape
        
        if size == 1:
            self.snapshots = data.copy()
            debug_print(f"Single processor mode. Data shape: {data.shape}")
            return
        
        else:
        
            debug_print(f"Data exchange started. Local shape: {data.shape}")

            # communicate the number of columns (snapshots) on each processor to everyone
            global_n_cols       = np.zeros(size, dtype=int) 
            global_n_cols[rank] = n_cols
            global_n_cols       = comm.allreduce(global_n_cols, op=MPI.SUM)
            
            debug_print(f"Global number of columns (snapshots): {global_n_cols}")
            
            # start sending data
            for target in range(size):
                
                # index of the vector segment to be sent
                i_s, i_e = distribute_mpi_work( n_rows, size, target )
                
                # extract the segment to be sent
                # slice must be copied to avoid issues with non-blocking send
                send_data = data[:,i_s:i_e].copy() 
                debug_print(f"Rank {rank} sending data to Rank {target}, shape: {send_data.shape}")
                req       = comm.Isend((send_data, MPI.FLOAT), dest=target, tag=target)
                debug_print(f"Rank {rank} sent data to Rank {target}")
                send_requests.append(req)
                
            
            # build the receive buffer
            
            i_s, i_e  = distribute_mpi_work( n_rows, size, rank )
            recv_buff = init_2Dflt_empty( np.sum(global_n_cols), i_e-i_s, 4)
            
            debug_print(f"Rank {rank} preparing receive buffer with shape: {recv_buff.shape}")
            
            # start receiving data
            for source in range(size):
                
                i_s, i_e = distribute_mpi_work( recv_buff.shape[0], size, source )
                debug_print(f"Rank {rank} receiving data from Rank {source}, shape: {recv_buff[i_s:i_e,:].shape}")
                req      = comm.Irecv((recv_buff[i_s:i_e,:], MPI.FLOAT), source=source, tag=rank)
                debug_print(f"Rank {rank} received data from Rank {source}, shape: {recv_buff[i_s:i_e,:].shape}")
                recv_requests.append(req)
            
            # wait for all sends and receives to complete
            
            debug_print("Waiting for all send requests to complete...")
            if send_requests:
                MPI.Request.Waitall(send_requests)
            debug_print("Waiting for all receive requests to complete...")
            if recv_requests:
                MPI.Request.Waitall(recv_requests)
                
            debug_print("All send and receive requests completed.")
            
            comm.barrier()
            
            self.snapshots = recv_buff
            
            debug_print(f"Data exchange complete. Received shape: {self.snapshots.shape}")
        
            return
        
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
#   - (Sayadi & Schmid, 2016)
# ----------------------------------------------------------------------

    def do_paradmd( self ):
        
        """
        No parameter needed.
        """
        
        # Check if the snapshots data are available
        
        self.N_t, self.len_snap_local = np.shape( self.snapshots )
  
        if self.N_t == 0 or self.N_t < 2 :
            
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
        
        Q1i, Ri = linalg.qr( Ai[:,:-1], mode='reduced' )
        
        
        # 2. Build the buffer also form the Rprime matrix
        
        Rprime =  init_1Dflt_empty( Ns*Ns*n_procs, self.kind )

            # ! be careful with the dimension of Ai(m,n) with m>n, otherwise,
            # you will get an error in broadcasting Ri into slice of Rprime.

        Rprime[ Ns*Ns*rank : Ns*Ns*(rank+1) ] = Ri.ravel()
        
        for r in range( n_procs ):
            
            self.comm.Bcast( Rprime[ Ns*Ns*r : Ns*Ns*(r+1) ], root=r )
        
        del Ri
        gc.collect()
        
        
        # 3. QR factorization of Rprime and get R, Q2i
        
        Q2i, R = linalg.qr( Rprime.reshape(Ns*n_procs, Ns), mode='reduced')

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
        
        Ur, S, Vt = linalg.svd( R, full_matrices=False )
        
        
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
        
        '''why C.transpose()?'''
        '''now transpose is not needed because of C order matrix'''
        
        B = np.array( np.matmul( C, Vt.conj().T ), order='F' )
        
        del C
        gc.collect()
        
        
        # 4. Compute S tilde
        
        Si_td = np.matmul( B, linalg.inv( sigma) )
        
        del B 
        gc.collect()
        
        
        # 5. Compute the eig(Si_st), with Si_d = Ui* x Ai(2=>Ns) x V x sigma^-1
        
        mu, Y = linalg.eig( Si_td )
        
        del Ai
        gc.collect()
        
#        print(mu)
        
        # ==============================================================
        # DMD modes in real space and amplitude
        # ==============================================================
        
        # DMD modes
        
        self.Phi_i = np.matmul( Ui, Y )
        
        
        # Build Vandermonde matrix
        Vand = init_2Dcmx_empty( Ns, Ns, self.kind*2 )
        
        Vand[:,0] = 1.0
        
        Vand[:,1] = mu

        for j in range(2, Ns): 
            Vand[:,j] = Vand[:,j-1] * mu
            
        
        # Build objective function matrices
        
        self.P = np.multiply( np.dot( Y.conj().T, Y )
                    , np.conj( np.dot( Vand, Vand.conj().T ) ) )
        
        self.q = np.conj( np.diag( np.dot( np.dot( np.dot( Vand, Vt.conj().T ), 
                                             sigma.conj().T ), Y ) ) )
        
        self.s = np.trace( np.dot( sigma.conj().T, sigma ) )

        
        # Find amplitudes directly
        
        alphas = linalg.solve( self.P, self.q )
        
        
        # Pass out variables
        
        self.mu = mu
        
        self.alphas = alphas
        
        self.freq = np.angle( mu )/self.dt/(2*np.pi)
        
        self.comm.barrier()
        
        if self.comm.rank == 0:
            print(f"DMD finished, got Phis shape {np.shape(self.Phi_i)}\n")

        sys.stdout.flush()


# ----------------------------------------------------------------------
# >>> Reconstruct from standard DMD                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/15  - created
#
# Desc
#   - reconstruct using all dmd modes
#   - root process collect and save the reconstructed data
#
# ----------------------------------------------------------------------

    def save_reconstruct( self ):
        
        """
        reconstruct using all dmd modes.\n
        root process collect and save the reconstructed data
        """
        
        # reconstruct on each process
        
        vand = np.vander( self.mu, 1, increasing=True ).astype(np.complex64)

        self.recons_data = self.Phi_i @ np.diag(self.alphas) @ vand
                
                
        # root process collect the reconstructed data of std dmd
        
        # - communicate between procs, same way as self.para_write_modes
        # - using len_mode_local ( from self.paradmd )
        
        self.len_mode = self.comm.allreduce( self.len_snap_local, op=MPI.SUM )
        
        send_counts = self.comm.gather( self.len_snap_local, root=0 )
        displace = [0] + list( np.cumsum( send_counts )[:-1] )
        
        
        if self.rank == 0:
            
            buf_recons = init_1Dcmx_empty( self.len_mode, self.kind*2 )
        
        else: buf_recons = None
        
        send_buf = self.recons_data.ravel().copy()
        
        self.comm.Gatherv( send_buf,
                           [buf_recons, send_counts, displace, MPI.COMPLEX8 ],
                           root=0 )
        
        # - save reconstructed data
        
        if self.rank == 0:
            
            recons_data_file = 'reconstructed_std_dmd.pkl'
            
            with open( recons_data_file, 'wb' ) as f:
                
                pickle.dump( buf_recons, f )
                
            del buf_recons
            gc.collect()
            
            print("reconstructed_std_dmd.pkl is saved.\n")
        
        
        
# ----------------------------------------------------------------------
# >>> Save P,q,s,Nt for further computing spdmd                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def save_Pqs( self ):
        
        """
        Root save Pqs.
        """
        
        # Check if P,q,s,N_t are available
        
        if self.P is None:
            
            raise ValueError("I don't have P,q,s now. Please do_paradmd() first"
                             " to compute objective functions.")
        
        
        # Save P,q,s,N_t with pickle
        
        with open( self.Pqs_file, "wb" ) as f:
            
            pickle.dump( self.P, f )
            
            pickle.dump( self.q, f )
            
            pickle.dump( self.s, f )
            
            pickle.dump( self.mu, f )
            
            pickle.dump( self.N_t, f )
            
            pickle.dump( self.freq, f )
            
            pickle.dump( self.alphas, f )
                

        print("Pqs.pkl file is written.\n")



# ----------------------------------------------------------------------
# >>> Read Pqs for spdmd                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_Pqs( self ):
        
        """
        Root read Pqs
        """
        
        # Check if Pqs.dat file is available
        
        if not ( os.path.exists( self.Pqs_file ) ):
            
            raise FileNotFoundError("Pqs.pkl is unavailable. Please run "
                                    "do_paradmd and save_Pqs first.")
   
        
        # Read P,q,s,N_t with pickle
        
        with open( self.Pqs_file,"rb") as f:
            
            self.P = pickle.load( f )
            
            self.q = pickle.load( f )
            
            self.s = pickle.load( f )
            
            self.mu = pickle.load( f )
            
            self.N_t = pickle.load( f )
            
            self.freq = pickle.load( f )
            
            self.alphas = pickle.load( f )
            
        
        print(f"{self.Pqs_file} is read.\n")



# ----------------------------------------------------------------------
# >>> Sparsity-promoting dmd                                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/11  - created
#
# Desc
#   - (Jovanović et al., 2014) 
#   - http://www.ece.umn.edu/users/mihailo//software/dmdsp/
#
# ----------------------------------------------------------------------

    def compute_spdmd( self, gamma=None, rho=None ):
        
        """
        gamma : default 10000. Larger gamma, less sp-mode.\n
        rho   : default 10000. Comparable to gamma, influence convergence speed.
        """
        
        # default parameters
        
        if gamma is None: gamma = 10000.
        
        if rho is None: rho = 10000.
        
        maxiter = 10000
        
        eps_abs = 0.000001
        
        eps_rel = 0.0001
        
        n = self.N_t - 1
        
        print(f"Matrix size is {n} x {n}.\n")
        
        # Selected mode index
        
        ind_spmode = None
        
        # Identity matrix
        
        I = np.eye( n )
        
        # Initialize data containers
        
        # N_none0    : number of non-zero amplitudes
        # Jsp        : square of F norm(before polish)
        # Jpol       : square of F norm(after polish)
        # Ploss      : optimal performance loss(after)
        # alphas_sp  : vector of amplitudes(before polish)
        # alphas_pol : vector of amplitudes(after polish)s
        
        N_none0 = 0     
        
        Jsp = 0.0          
        
        Jpol = 0.0         
        
        Ploss = 0.0        
        
        alphas_sp = np.zeros( n, dtype=complex )  
        
        alphas_pol = np.zeros( n, dtype=complex ) 
    
        
        # Check if P, q, s are available
        
        if self.P is None: raise ValueError('Please compute dmd first.')

        
        # Cholesky factorization of matrix P + (rho/2)*I
        
        Prho = ( self.P + (rho/2.0) * np.eye( n ) )
        
        Plow = linalg.cholesky( Prho )
        
        
        # Initial condition
        
        y = np.zeros( n )
        z = np.zeros( n )
        
        # Use ADMM to solve the gamma-parameterized problem
        
        for iter in range( maxiter ):
            
            t_ref = time.time()
        
            # x-minimization step 
            # x =  (P + (rho/2)*I)\(q + (rho)*u)
            
            u = z - ( 1.0/rho ) * y 
            xnew = linalg.solve( Plow.conj().T,
                        linalg.solve( Plow, (self.q + (rho/2.0)*u )))
        
        
            # z-minimization step
            
            a = ( gamma/rho ) * np.ones( n )
            v = xnew + ( 1.0/rho ) * y
            
            # Soft-thresholding of v
            
            znew = ( ( 1.0 - a/np.abs(v) ) * v ) * (np.abs(v) > a)
            
            
            # Primal and dual residuals
            
            res_prim = linalg.norm( xnew - znew )
            res_dual = rho * linalg.norm( znew - z )
            
            
            # Lagrange multiplier update step
            
            y = y + rho * ( xnew - znew )
            
            
            # Stopping criteria
            
            eps_prim = np.sqrt( n ) * eps_abs + \
                        eps_rel * max( linalg.norm(xnew), linalg.norm(znew) )

            eps_dual = np.sqrt( n ) * eps_abs + eps_rel * linalg.norm( y )
            
            print( f"Iter. {iter+1:5d} took {time.time()-t_ref:5.2f}s.",end='')
            print( f" - res_prim : {res_prim:10.2e} < {eps_prim:10.2e}",end='') 
            print( f" - res_dual : {res_dual:10.2e} < {eps_dual:10.2e}",end='') 
            print( f" Current nr. of selected modes:{np.count_nonzero(z):5d}")
            
            
            if (res_prim < eps_prim) and (res_dual < eps_dual): break
        
            else: z = znew
        
        # verbose after ADMM iteration
        
        print( "\nIteration has finished, solving KKT system ... \n" )
        
        
        # Record output data 
        
        alphas_sp[:] = z
        
        N_none0 = np.count_nonzero( alphas_sp[:] )
        
        Jsp = np.real( z.conj().T @ self.P @ z ) \
                    - 2.0 * np.real( self.q.conj().T @ z) + self.s
                    
        
        # Polishing of the nonzero amplitudes
        # Form the constraint matrix E for E^T x = 0
        
        ind_zero = np.where( np.abs(z) < 1.e-12 )[0]
        
        ind_spmode = np.where( np.abs(z) > 1.e-12 )[0] 
        
        m = len( ind_zero )
        
        E = I[:,ind_zero]
        
        
        # Form KKT system for the optimality conditions
        
        KKT = np.block( [[self.P, E], [E.T, np.zeros((m,m))]] )
        
        rhs = np.block( [self.q, np.zeros(m)] )
        
        
        # Solve KKT system 
        
        sol = linalg.solve( KKT, rhs )
        
        
        # Vector of polished (optimal) amplitudes
        
        alphas_pol[:] = sol[0:n]
        
        # Polished (optimal) least-squares residual
        
        Jpol = np.real( sol[0:n].conj().T @ self.P @ sol[0:n] )\
                    - 2.0 * np.real(self.q.conj().T @ sol[0:n]) + self.s

        # Polished (optimal) performance loss
        
        Ploss = 100.0 * np.sqrt( Jpol/self.s )
            
        
        # Record output data
        
        self.Jsp = Jsp
        self.Jpol = Jpol
        self.Ploss = Ploss
        self.alphas_sp = alphas_sp
        self.alphas_pol = alphas_pol
        self.ind_spmode = ind_spmode

        # Compute other parameters
        
        # - beta : growth rate
        # - psi  : relative amplitudes
        
        self.beta = np.log( np.abs( self.mu ) )
        
        self.psi_pol = np.abs( self.alphas_pol )/max( np.abs(self.alphas_pol) )
        
        self.psi = np.abs( self.alphas )/max( np.abs(self.alphas_pol) )
       


# ----------------------------------------------------------------------
# >>> Save result of spdmd for 2nd round DMD                     (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/04  - created
#
# Desc
#   - save ind_spmode, alphas_sp, alphas_pol
#
# ----------------------------------------------------------------------

    def save_spdmd_result( self ):
        
        """
        No parameter needed. \n
        Save spdmd result into self.spdmd_result_file
        """
        
        filename = self.spdmd_result_file
        
        if hasattr(self, 'ind_spmode'):
        
            with open( filename, "wb") as f:
                
                pickle.dump( self.ind_spmode, f )
                
                pickle.dump( self.alphas, f )
                
                pickle.dump( self.alphas_sp, f )
                
                pickle.dump( self.alphas_pol, f )
                
                pickle.dump( self.Ploss, f )
                
                pickle.dump( self.Jsp, f )
                
                pickle.dump( self.Jpol, f )
                
                pickle.dump( self.beta, f )
                
                pickle.dump( self.psi_pol, f )
                
                pickle.dump( self.psi, f )
                
                pickle.dump( self.St, f )
                
                pickle.dump( self.mu, f )
                
        else: raise ValueError("Please compute spdmd first!")

        print("spdmd results are saved.\n")

# ----------------------------------------------------------------------
# >>> Save index of sparsity-promoting modes                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/03  - created
#
# Desc
#   - save the index of sparsity-promoting modes to csv.
#   - Just for manually check, will not be used.
# ----------------------------------------------------------------------

    def save_ind_spmode( self ):
        
        """
        No parameter needed. \n 
        Save the index of sparsity-promoting modes to ind_spmode.csv
        """
        
        if hasattr(self, 'ind_spmode'):
        
            df = pd.DataFrame( self.ind_spmode, columns=['ind_spmode'])
            
            df['alphas'] = self.alphas[ tuple([self.ind_spmode]) ]
            
            df['alphas_pol'] = self.alphas_pol[ tuple([self.ind_spmode]) ]
            
            df['alphas_sp'] = self.alphas_sp[ tuple([self.ind_spmode]) ]
            
            df['beta'] = self.beta[ tuple([self.ind_spmode]) ]
            
            df['freq'] = self.freq[ tuple([self.ind_spmode]) ]
            
            df['psi'] = self.psi[ tuple([self.ind_spmode]) ]
            
            df['psi_pol'] = self.psi_pol[ tuple([self.ind_spmode]) ]
            
            if self.St is not None:
                
                df['St'] = self.St[ tuple([self.ind_spmode]) ]
            
            else: 
                
                print("\nWarning: St is not computed when saving ind_spmode\n")
            
            # make a copy of sorted( based on St ) dataframe
            
            df_unsorted = df.copy()
            
            # align each column
            
            df = df.drop( df[ df['St']<0.0 ].index )
            df = df.sort_values(by='St')
            
            df = df.astype(str)
            
            max_lengths = df.map(len).max()
            
            for column in df.columns:
                df[column] = df[column].apply(lambda x:
                                              x.ljust(max_lengths[column]))
            
            df.to_csv( 'ind_spmode.csv', 
                        sep='\t', 
                        index=False )
            
            # output the unsorted, full index
            df_unsorted = df_unsorted.astype(str)
            
            max_lengths = df_unsorted.map(len).max()
            
            for column in df_unsorted.columns:
                df_unsorted[column] = df_unsorted[column].apply(lambda x:
                                              x.ljust(max_lengths[column]))
            
            # output to csv
            
            df_unsorted.to_csv( 'ind_spmode.csv', 
                                sep='\t', 
                                index=False,
                                mode='a')
            
        
        else: raise ValueError("Please compute spdmd first!")
        
        print("Indexes of selected sparsity promoting modes are saved.\n")



# ----------------------------------------------------------------------
# >>> Read spdmd_result file and broadcast                          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/03  - created
#
# Desc
#
# ----------------------------------------------------------------------
    
    def read_spdmd_result( self ):
        
        """
        No parameter needed. \n
        Read spdmd_result from self.spdmd_result_file
        """
        
        self.ind_spmode = None
        self.alphas_sp  = None
        self.alphas_pol = None
        
        if self.rank == 0:
            
            filename = self.spdmd_result_file
            
            with open( filename, "rb" ) as f:
            
                self.ind_spmode = pickle.load( f )
                
                self.alphas     = pickle.load( f )
                
                self.alphas_sp  = pickle.load( f )
                
                self.alphas_pol = pickle.load( f )
                
                self.Ploss      = pickle.load( f )
                
                self.Jsp        = pickle.load( f )
                
                self.Jpol       = pickle.load( f )
                
                self.beta       = pickle.load( f )
                
                self.psi_pol    = pickle.load( f )
                
                self.psi        = pickle.load( f )
                
                self.St         = pickle.load( f )
                
                self.mu         = pickle.load( f )
                
            print(col.fg.green,"Finish reading in spdmd results.",
                  col.reset,end='\n')
            
            print(f"\nIndexes of selected modes are:\n{self.ind_spmode}")
            
            print(f"Ploss = {self.Ploss:8.3f} % .  ",end='')
            
            print(f"Jsp = {self.Jsp:10.3e}, Jpol = {self.Jpol:10.3e}.\n")
        
        
        self.ind_spmode = self.comm.bcast( self.ind_spmode, root=0 )
        
        self.alphas_sp = self.comm.bcast( self.alphas_sp, root=0 )
        
        self.alphas_pol = self.comm.bcast( self.alphas_pol, root=0 )



# ----------------------------------------------------------------------
# >>> parallel write selected dmd modes                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/03  - created
#
# Desc
#
# ----------------------------------------------------------------------
    
    def para_write_modes( self ):
        
        """
        No parameter needed. \n
        Parallel write selected dmd modes into mode_xxxxx.pkl
        """
        
        # Enter the existing directory or make a new one
        
        if os.path.exists( self.dmdmodes_dir ):
            
            os.chdir( self.dmdmodes_dir )
        
        else:
            
            if self.rank == 0:
                
                os.mkdir( self.dmdmodes_dir )
                
                print("dmdmodes dir does NOT exist, new dir is made.\n")
        
            self.comm.barrier()
            
            os.chdir( self.dmdmodes_dir )
        
        
        # Collect modes and write by root
        
        # len_mode is the length of a whole flow field,
        # on contrary, len_snap_local is the length of snapshot on each proc
        
        self.len_mode = self.comm.allreduce( self.len_snap_local, op=MPI.SUM )
        
        
        # Gather mode one by one from other procs to root
        # Parameters for MPI.Gatherv: 
        # send_counts: a list of number of data elements on each proc
        # displace: a list of locations where gathered data shoud put
        
        send_counts = self.comm.gather( self.len_snap_local, root=0 )
        displace = [0] + list( np.cumsum( send_counts)[:-1] )

        for i in range( len( self.ind_spmode ) ):
         
#            print(f"send_counts is {send_counts}, rank is {self.rank}")           
#            print(f"displace is {displace},rank is {self.rank}")
            
            if self.rank == 0:
                
                buf_mode = init_1Dcmx_empty( self.len_mode, self.kind*2 )
            
            else: buf_mode = None
                
            indx = self.ind_spmode[i]
            
            # Phi_i is a tall matrix,
            # deep copy to make data contagious for mpi
            
            send_buf = self.Phi_i[:,indx].copy()
                       
            self.comm.Gatherv( send_buf, 
                               [buf_mode, send_counts, displace, MPI.COMPLEX8 ],
                               root=0 )
            
            # Save modes
            
            if self.rank == 0:
                
                mode_filename = "mode_" + f"{indx:05d}" + ".pkl"
                
                with open( mode_filename, 'wb' ) as f:
                    
                    # St is not necessary actually, can be derived from mu.
                    
                    print(f'mode {i}')
                    print(f"indx is {indx}")
                    print(f"alpha is {self.alphas[indx]}")
                    print(f"alpha_sp is {self.alphas_sp[indx]}")
                    print(f"alpha_pol is {self.alphas_pol[indx]}")
                    print(f"mu is {self.mu[indx]}")
                    print(f"St is {self.St[indx]}")
                    print("=============\n")
                    
                    pickle.dump( indx, f )
                    pickle.dump( self.alphas[indx], f )
                    pickle.dump( self.alphas_sp[indx], f )
                    pickle.dump( self.alphas_pol[indx], f )
                    pickle.dump( self.mu[indx], f )
                    pickle.dump( self.St[indx], f )
                    pickle.dump( buf_mode, f )

                    del buf_mode
                    del send_buf
                    gc.collect()



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
