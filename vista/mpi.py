# -*- coding: utf-8 -*-
'''
@File    :   mpi.py
@Time    :   2024/10/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import sys
from   mpi4py            import MPI


class MPIenv:

    def __init__(self):
        self._comm = MPI.COMM_WORLD
        self._rank = self.comm.Get_rank()
        self._size = self.comm.Get_size()
        self._root = 0
        
    @property
    def rank(self): return self._rank

    @property
    def comm(self): return self._comm

    @property
    def size(self): return self._size

    @property
    def root(self): return self._root

    @property
    def is_root(self): return self._rank == self._root

    def barrier(self): self._comm.Barrier()

# ----------------------------------------------------------------------
# >>> master distribute                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def master_distribute(self, tasks):
        
        """
        tasks: list of tasks to be distributed to workers
        """
        
        # check the allocated 
        
        if self.is_root:
            
            task = tasks
            task_index = 0

            num_workers = self.size - 1
            
            if num_workers == 0:
                print("No workers available. Master should do all the tasks.")
                return
            
            # number of closed workers
            closed_workers = 0   

            while closed_workers < num_workers:
                
                status = MPI.Status()
                
                # receive a request(None) from any worker, and get the source tag
                
                data   = self.comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
                source = status.Get_source()  
                
                # if there are still tasks to be distributed
                # send the task to the worker who just contacted the master
                
                if task_index < len(tasks): 
                    task = tasks[task_index]
                    print(f"Master sending task {task} to worker {source}.  {task_index}/{len(tasks)}.")
                    self.comm.send(task_index, dest=source, tag=0)
                    task_index += 1
                
                # if no task left, send termination signal to the worker
                
                else:
                    print(f"Master sending termination signal to worker {source}")
                    self.comm.send(None, dest=source, tag=1) 
                    closed_workers += 1
            
                sys.stdout.flush()


# ----------------------------------------------------------------------
# >>> worker receive                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def worker_receive( self ):
        
        """
        receive task from master.
        
        -return: task, if None, then please exit.
        """
        
        # send request to master
        
        self.comm.send(None, dest=0)
        
        # receive task from master
        
        task_index = self.comm.recv(source=0, tag=MPI.ANY_TAG)
        
        # if task is None, then exit
        
        if task_index is None:
            print(f"Worker {self.rank} received termination signal. Exiting.")
        
        else:
            print(f"Worker {self.rank} is processing task {task_index}")
        
        sys.stdout.flush()
        
        return task_index


# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    mpi = MPIenv()

    tasks = ['task1', 'task2','task3','task4','task5']
    
    def do_task( task ):
        print(f"rank_{mpi.rank} is doing {task}.")

    if mpi.size == 1:
        print("No workers available. Master should do all the tasks.")
    
        for task in tasks:
            do_task( task )
    
    else: 
        if mpi.rank == 0:
            mpi.master_distribute( tasks )

        else:
            while True:
                task_index = mpi.worker_receive()
                
                if task_index is None: break
                else: do_task( tasks[task_index] )


# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/10/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()

