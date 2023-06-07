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

import numpy             as     np

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
# - Only add files when key is matched
# ----------------------------------------------------------------------

def get_filelist( FoldPath, key=None ):
    
    FileList = []
    
    for home, dirs, files in os.walk( FoldPath ):
        
        for filename in files:
            
            # Filelist need to contain the whole path
            if key is not None:
                
                if key in filename:
                    
                    FileList.append( os.path.join(home,filename) )
            else:
                
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
#   5 : 1221 case, D/delta = 0.125
#
# - wave length lambda = 1.04
#
# ----------------------------------------------------------------------

def is_above_wavywall( y, z, Case = 1):
    
    A = 0.26
    
    len_w = 1.04
    
    if   Case == 1:   D = 10.4
    elif Case == 2:   D = 5.2    
    elif Case == 3:   D = 2.6       
    elif Case == 4:   D = 1.3
    # Case 5: the wavelength is smaller and do not have flat valley
    elif Case == 5:
        len_w = 0.65
        D     = 0.65
    
    z0 = z % D 
    
    # when z0 is less than half wave length
    if abs(z0) <= (len_w*0.5):
    
        y_w = -A + A*math.cos( z0/len_w*2*math.pi )
    
    elif abs(z0-D) <= (len_w*0.5):
    
        y_w = -A + A*math.cos( (z0-D)/len_w*2*math.pi )
    
    else:
        # for case 5, this will not happen
        y_w = -2*A
    
    # compare wall with cell center location
    if y >= y_w: 
    
        above_wall = True
    
    else:
    
        above_wall = False
    
    return above_wall

# ----------------------------------------------------------------------
# >>> Mean of a list                                            ( 4 )
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
# - simple function applied in pandas 
#
# ----------------------------------------------------------------------

def mean_of_list(lst):
    return sum(lst) / len(lst)


# ----------------------------------------------------------------------
# >>> To dictionary                                             ( 5 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/18  - created
#
# Desc
#   - match two list(vectors) to a dictionary
#
# ----------------------------------------------------------------------

def to_dictionary( keys, values ):
    return { key:value for key, value in zip(keys, values) }



# ----------------------------------------------------------------------
# >>> Read parameter from file                                   (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def read_case_parameter( filename ):
    
    # create a dictionary for parameter pairs
    
    parameters = {}


    # read file content line by line
    
    with open( filename, 'r') as f:
        
        lines = f.readlines()

        for line in lines:
            
            # ignore notation line and space line
            
            if line.strip() and not line.startswith("#"):
                
                key, value = line.strip().split("=")
                
                # remove extra space and comma
                
                key = key.strip()
                value = value.strip().rstrip(",")
                
                # add parameter pair into the dictionary
                
                parameters[key] = value
                
    
    return parameters
                
                
                
# ----------------------------------------------------------------------
# >>> Main: for testing and debugging                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/15  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":
    
    filename = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/snapshots/snapshot_test_z/case_parameters"
    
    parameters = read_case_parameter(filename)
    
    for key, value in parameters.items():
        print(f"key is {key}, value is {value}")