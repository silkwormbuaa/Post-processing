#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_tool.py
@Time    :   2022/10/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Some general tools 
'''

import os 


# ----------------------------------------------------------------------
# >>> GET FILE LIST                                               ( 0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/13  - created
#
# Desc
# 
# - All files, including files in any subfolder will be retrieved
#
# ----------------------------------------------------------------------

def get_filelist( FoldPath ):
    
    FileList = []
    
    for home, dirs, files in os.walk( FoldPath ):
        
        for filename in files:
            
            # Filelist need to contain the whole path
            
            FileList.append( os.path.join(home,filename) )
            
    return FileList

