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

class Logger:
    
    def __init__( self, filename=None ):
        
        if filename is None:
            filename = "logfile.log"
            
        self.terminal = sys.stdout
        name, _ = os.path.splitext( filename )
        self.log = open(f"{name}.log", "a")
   
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass    

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