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


# ----------------------------------------------------------------------
# >>> IF_OVERLAP                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/15  - created
#
# Desc
#
# - check if two rectangular region overlap with each other
#
# Input
# 
# - two lists which contain (xmin,ymin,xmax,ymax) respectively
#
# ----------------------------------------------------------------------

def if_overlap( rect1, rect2 ):    
    
    notOverlap = ( rect1[2] <= rect2[0] ) or \
                 ( rect1[0] >= rect2[2] ) or \
                 ( rect1[3] <= rect2[1] ) or \
                 ( rect1[1] >= rect2[3] )
    
    Overlap = not notOverlap
    
    return Overlap


