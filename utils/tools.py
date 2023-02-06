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

import math


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
    
    FileList.sort()
            
    return FileList


# ----------------------------------------------------------------------
# >>> IF_OVERLAP                                                ( 1 )
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


# ----------------------------------------------------------------------
# >>> If Segment Penetrate Zone                                  ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/24  - created
#
# Desc
#
# - Check if a streamwise segment penetrate a zone
#
# ----------------------------------------------------------------------

def if_penetrate( zone_range, segment ):
    
    Penetrate = ( segment[0] < zone_range[3] ) and \
                ( segment[0] > zone_range[1] ) and \
                ( segment[1] < zone_range[0] ) and \
                ( segment[2] > zone_range[2] )
    
    return Penetrate

# ----------------------------------------------------------------------
# >>> If a point is above wall                                   ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/06  - created
#
# Desc
#
# - special routine for wavy wall case to check if a point is 
#   above or below wall
#
# - select which wavy wall case
#
#   1 : 1014 case, D/delta = 2
#   2 : 0926 case, D/delta = 1
#   3 : 0825 case, D/delta = 0.5
#   4 : 0927 case, D/delta = 0.25
#
# - wave length lambda = 1.04
#
# ----------------------------------------------------------------------

def is_above_wavywall( z, y, case = 1):
    
    A = 0.26
    
    len_w = 1.04
    
    if   case == 1:   D = 10.4
    elif case == 2:   D = 5.2    
    elif case == 3:   D = 2.6       
    elif case == 4:   D = 1.3
    
    z0 = z % D 
    
    # when z0 is less than half wave length
    if abs(z0) <= (len_w*0.5):
    
        y_w = -A + A*math.cos( z0/len_w*2*math.pi )
    
    elif abs(z0-D) <= (len_w*0.5):
    
        y_w = -A + A*math.cos( (z0-D)/len_w*2*math.pi )
    
    else:
    
        y_w = -2*A
    
    # compare wall with cell center location
    if y >= y_w: 
    
        above_wall = True
    
    else:
    
        above_wall = False
    
    return above_wall


