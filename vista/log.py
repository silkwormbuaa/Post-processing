# -*- coding: utf-8 -*-
'''
@File    :   log.py
@Time    :   2023/10/27 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import time
import atexit

class Logger:
    
    def __init__( self, filename=None ):
        
        if filename is None:
            filename = "logfile.log"
            
        self.terminal = sys.stdout
        name, _ = os.path.splitext( filename )
        self.log = open( os.path.join( os.getcwd(), f"{name}.log"), "a")
        
        sys.stdout = self
        atexit.register(self._on_exit)
   
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass

    def _on_exit(self):
        
        # write time in year-month-day hour:minute:second format
        exit_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        self.log.write(f"\n=== Program exited at {exit_time} ===\n\n")
        self.log.close()

# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/27  - created
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
# 2023/10/27  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()