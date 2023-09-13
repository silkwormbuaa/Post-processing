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

import shutil

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
# >>> copy_file_with_same_subfolder                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/26  - created
#
# Desc
#   - a function to copy file with same subfolder
#     like copying a kind of snapshot under snapshots/snapshot_xxx/xxx 
# ----------------------------------------------------------------------

def copy_file_with_same_subfolder( src_folder, dest_folder, key ):
    
    src_file_list = get_filelist( src_folder, key)

    dest_file_list = []

    for file in src_file_list:
        
        rear = file.split(src_folder,maxsplit=1)[1]
        
        dest_file = dest_folder.rstrip('/')+'/'+rear.lstrip('/')
        
        dest_file_list.append(dest_file)

    n_files = len(dest_file_list)

    for i in range(n_files):
        
        if os.path.exists(dest_file_list[i].split(key)[0]):
        
            shutil.copy2(src_file_list[i],dest_file_list[i])
            
            if i%(n_files//100) == 0:
                
                print(f"copy progress : {int(i/n_files*100)}%")
        
        else: raise FileNotFoundError(dest_file_list[i].split(key)[0])


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

def if_overlap_2d( rect1, rect2 ):    
    
    notOverlap = ( rect1[2] <= rect2[0] ) or \
                 ( rect1[0] >= rect2[2] ) or \
                 ( rect1[3] <= rect2[1] ) or \
                 ( rect1[1] >= rect2[3] )
    
    Overlap = not notOverlap
    
    return Overlap


# ----------------------------------------------------------------------
# >>> IF OVERLAP 3D                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/01  - created
#
# Desc
#
#   - check if two box overlap
#   - bbox1 = [x1min,x1max,y1min,y1max,z1min,z1max]
#   - bbox2 = [x2min,x2max,y2min,y2max,z2min,z2max]
# ----------------------------------------------------------------------

def if_overlap_3d( bbox1, bbox2 ):    
    
    notOverlap = ( bbox1[1] <= bbox2[0] ) or \
                 ( bbox1[0] >= bbox2[1] ) or \
                 ( bbox1[3] <= bbox2[2] ) or \
                 ( bbox1[2] >= bbox2[3] ) or \
                 ( bbox1[5] <= bbox2[4] ) or \
                 ( bbox1[4] >= bbox2[5] )
    
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
#
# - wave length lambda = 1.04
#
# ----------------------------------------------------------------------

def is_above_wavywall( y, z, Case = 1):
    
    """
    y: float
    z: float
    case=1 : 1014 case, D/delta = 2
    case=2 : 0926 case, D/delta = 1
    case=3 : 0825 case, D/delta = 0.5
    case=4 : 0927 case, D/delta = 0.25
    case=5 : 1221 case, D/delta = 0.125
    """
    
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
# >>> define wall shape                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

def define_wall_shape( z_list:np.array, Case = 1, write=True, scale=True ):

    """
    z_list: np.array of z coordinates
    case=1 : 1014 case, D/delta = 2
    case=2 : 0926 case, D/delta = 1
    case=3 : 0825 case, D/delta = 0.5
    case=4 : 0927 case, D/delta = 0.25
    case=5 : 1221 case, D/delta = 0.125
    
    write out wall_X.dat or return list of y
    """
       
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
    
    y_list = []
    
    for z in z_list:
    
        z0 = z % D 
        
        # when z0 is less than half wave length
        if abs(z0) <= (len_w*0.5):
        
            y_w = -A + A*math.cos( z0/len_w*2*math.pi )
        
        elif abs(z0-D) <= (len_w*0.5):
        
            y_w = -A + A*math.cos( (z0-D)/len_w*2*math.pi )
        
        else:
            # for case 5, this will not happen
            y_w = -2.0*A

        y_list.append( y_w )
    
    if write:
        
        zycor = np.column_stack((z_list,y_list))
        
        if scale: zycor = zycor/5.2
        
        np.savetxt( "wall_X.dat",
                    zycor,
                    fmt="%15.7f",
                    delimiter=" ",
                    comments="",
                    header='x    y')
    
    else: 
        
        return np.array( y_list )
    
    

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

    if not os.path.exists( filename ):
        raise FileExistsError("Please set case_parameters file!")

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
# >>> Distribute MPI work                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/04  - created
#
# Desc
#
#   - Distribute the tasks each processor should take 
#   - return values are 'i_start' and 'i_end' representing the start and 
#     end index of a list of N tasks starting from 0. (0,1,...,N-1)
#   - if there is no task, i_start = None
# ----------------------------------------------------------------------

def distribute_mpi_work(n_tasks, n_processors, rank):
    
    # algorithm divides number of blocks in a dealing way
    
    n_tasks_local = n_tasks // n_processors       # floor dividing
    
    left_tasks = n_tasks - n_tasks_local*n_processors
    
    
    if rank < left_tasks:
        
        n_tasks_local = n_tasks_local + 1
        
        i_start = 0 + rank*n_tasks_local
        
        i_end = i_start + n_tasks_local
    
    
    elif n_tasks_local > 0:
        
        i_start = 0 + rank*n_tasks_local + left_tasks
        
        i_end = i_start + n_tasks_local    
    
        
    else:
        
        n_tasks_local = 0
        
        i_start = None
        
        i_end = None
        
    
    return i_start, i_end
        
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