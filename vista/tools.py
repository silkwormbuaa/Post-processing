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
import sys
import cv2
import math
import shutil
import numpy             as     np
import matplotlib.pyplot as     plt
from   scipy.interpolate import interp1d
from   scipy.optimize    import root

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
    
    # check if FoldPath exists
    if os.path.exists( FoldPath ):
        
        for home, dirs, files in os.walk( FoldPath ):
            
            for filename in files:
                
                # Filelist need to contain the whole path
                if key is not None:
                    
                    if key in filename:
                        
                        FileList.append( os.path.join(home,filename) )
                else:
                    
                    FileList.append( os.path.join(home,filename) )
        
        FileList.sort()
    
    else:
        print(f"Folder {FoldPath} does not exist!"); sys.exit()
            
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
# >>> point_in_box                                          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

def point_in_box( xyz, box ):
    
    """
    xyz: [x,y,z] of the given point
    box: [xmin,ymin,zmin,xmax,ymax,zmax] of the 3D box
    """
    
    if ( xyz[0] >= box[0] ) and ( xyz[0] <= box[3] ) and \
       ( xyz[1] >= box[1] ) and ( xyz[1] <= box[4] ) and \
       ( xyz[2] >= box[2] ) and ( xyz[2] <= box[5] ):
        
        return True
    
    else:
        return False
    
    
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

def define_wall_shape( z_list:np.array, Case = None, casecode=None,
                       yshift=None, write=True, scale=True ):

    """
    z_list: np.array of z coordinates
    case=1 : 221014 case, D/delta = 2.0,  A = 0.26
    case=2 : 220926 case, D/delta = 1.0,  A = 0.26
    case=3 : 220825 case, D/delta = 0.5,  A = 0.26
    case=4 : 220927 case, D/delta = 0.25, A = 0.26
             231124 case, D/delta = 0.25, A = 0.26, Re_tau = 1226
    case=5 : 221221 case, D/delta = 0.125,A = 0.26
    case=6 : 240211 case, D/delta = 0.25, A = 0.13
    case=7 : 240210 case, D/delta = 0.25, A = 0.52
             241018 case, D/delta = 0.25, A = 0.52, Re_tau = 1226
    case=8 : 241017 case, D/delta = 0.25, A = 0.06837, Re_tau = 1226
    case=9 : 241030 case, D/delta = 0.0625, A = 0.068395, Re_tau = 1226
    
    write out wall_X.dat or return list of y
    """
       
    A = 0.26
    
    len_w = 1.04
    
    if casecode is not None:
        if   casecode == '221014': Case=1
        elif casecode == '220926': Case=2
        elif casecode == '220825': Case=3
        elif casecode == '220927': Case=4
        elif casecode == '221221': Case=5
        elif casecode == '240211': Case=6
        elif casecode == '240210': Case=7
        elif 'smooth' in casecode: Case=0
        elif casecode == '231124': Case=4
        elif casecode == '241017': Case=8
        elif casecode == '241018': Case=7
        elif casecode == '241030': Case=9
        elif casecode == '250120': Case=2
        elif casecode == '250218': Case=2
        elif casecode == '250304': Case=2
        
    
    if   Case == 1:   D = 10.4
    elif Case == 2:   D = 5.2    
    elif Case == 3:   D = 2.6       
    elif Case == 4:   D = 1.3
    # Case 5: the wavelength is smaller and do not have flat valley
    elif Case == 5:   len_w = 0.65; D = 0.65
    elif Case == 6:   D = 1.3;      A = 0.13
    elif Case == 7:   D = 1.3;      A = 0.52
    elif Case == 8:   D = 1.3;      A = 0.06837
    elif Case == 9:   D = 0.325;    A = 0.068395;  len_w = 0.27359
    
    y_list = []
    
# - define wall shape    

    if Case == 0:
        y_list = np.zeros_like(z_list)
    
    else:
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

# - shift y?
 
    if yshift is not None:
        y_list = np.array(y_list) + yshift
        
    
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
# >>> bilinear interpolation                                     (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/22  - created
#
# Desc
#
#      f(3)  # . . . . . . . . # f(4)
#     (x1,y2).                 . (x2,y2)
#            .                 .
#            .         f(x,y)  .
#            .            *    .
#            .                 .
#      f(1)  # . . . . . . . . # f(2)
#       (x1,y1)                (x2,y1)
#
# ----------------------------------------------------------------------

def bilin_interp(x1,x2,y1,y2,f,x,y):
    
    
    """
    x1,x2,y1,y3 : coordinates at left-bottom and right-top corner
    f           : values at (x1,y1)--(x2,y1)--(x1,y2)--(x2,y2)
    x,y         : interpolation point location 
    
    return : f(x,y)
    """

#---- calculate weights

    u = (x - x1) / (x2 - x1)
    v = (y - y1) / (y2 - y1)

#---- calculate interpolated value

    value =   (1.0-u)*(1.0-v) * f[0] \
            +      u *(1.0-v) * f[1] \
            + (1.0-u)*     v  * f[2] \
            +      u *     v  * f[3]
    
    return value


# ----------------------------------------------------------------------
# >>> find indices                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

def find_indices( arr:np.ndarray, target:float, mode='bisection' ):
    
    """
    arr    : 1D list or numpy array, must monotically increase
    target : target value
    mode   : 'bisection' or 'sequential'  
    
    return : left and right indices 
    """

 
    if mode == 'bisection':
        
        # check if the target is out of bound    
        if target < arr[0] or target > arr[-1]:
            raise ValueError("The target value is out of bound.")
    
        left, right = 0, len(arr) - 1

        while left <= right:
            mid = left + (right - left) // 2  # avoid integer overflow

            if arr[mid] == target:
                return mid, mid + 1  # find target value
            elif arr[mid] < target:
                left = mid + 1
            else:
                right = mid - 1

        # return left and right index
        return right, left

    if mode == 'sequential':
        
        left = 0
        
        while arr[left] < target:
            
            left += 1
        
        return left -1 , left


# ----------------------------------------------------------------------
# >>> convert image format                                       (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/18  - created
#
# Desc
#
# ----------------------------------------------------------------------

def convert_image_format(input_file, 
                         output_file, 
                         output_format, 
                         outputsize=None):
    """
    input_file  :  the input image
    output_file :  the output image
    output_format : 'jpg','pdf','eps','png','bmp'
    outputsize : (width,height)
    
    ! But actually now can only convert from bitmap to vector image, cannot the
    other way around.
    
    ! Suggest generate all format figures at the beginning. When there is really 
      need to convert figures, use linux command 'convert'
    """
    try:
        
        # open image
        img = plt.imread( input_file )
            
        # resize image
        if outputsize:
            img = img.resize(outputsize)
        
        # get output path
        output_dir, output_filename = os.path.split(output_file)
        output_name, _ = os.path.splitext(output_filename)
        output_file = os.path.join(output_dir, f"{output_name}.{output_format}")

        plt.imshow(img)
        plt.axis('off')
        
        # convert and save figure
        plt.savefig( output_file, dpi=1000 )
        print(f"{input_file} has been converted as {output_file}")

    except Exception as e:
        print(f"Convert failed:{e}")
        

# ----------------------------------------------------------------------
# >>> Linear-growing grid vector                                 (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/19  - created
#
# Desc
#
# ----------------------------------------------------------------------

def lin_grow(xs,ds,expratio,len=None,upbound=None):
    
    """
    xs       : starting point of array \n
    ds       : starting spacing of array \n
    expratio : expansion ratio of spacing \n
    len      : total length of cells ( nr. of grid points - 1) \n
    upbound  : upper bound of grid point \n
    
    return:    an array of linear-growing spacing grid x and dx array
    """
    # if give the total length of array
    
    if len is not None:
        x = np.zeros(len+1,dtype='f8')
        x[0] = xs
        x[1] = xs+ds
        
        for i in range(2,len+1):
            ds = ds*expratio
            x[i] = x[i-1] + ds
        
        d = np.diff(x)
    
    
    # if give the upper bound of array
    
    elif upbound is not None:
        x = []
        x.append(xs)
        x.append(xs+ds)
        
        while x[-1] < upbound:
            ds = ds*expratio
            x.append(x[-1]+ds)
            
        x[-1] = upbound
        
        d = np.diff(x)
    
    else: raise ValueError("Please give 'len' OR 'upbound' target array.")
    
    return x, d


# ----------------------------------------------------------------------
# >>> interpolator1d                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/06/20  - created
#
# Desc
#
# ----------------------------------------------------------------------

def create_linear_interpolator(x, y):
    """
    create a linear interpolator function based on the given x and y data.

    inputs:
    x : array-like
    y : array-like

    return:
    function
    """
    # create a linear interpolator
    linear_interp = interp1d(x, y, kind='linear', fill_value="extrapolate")

    def interpolator(x0):
        """
        根据给定的 x0 值进行线性插值，返回对应的 y0 值。

        参数:
        x0 : float or array-like
            需要插值的 x 坐标值。

        返回:
        float or array-like
            对应的插值后的 y 坐标值。
        """
        return linear_interp(x0)

    return interpolator


# ----------------------------------------------------------------------
# >>> crop image (RGB image only)                                            
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/09/25  - created
#
# Desc
#
# ----------------------------------------------------------------------

def crop_border(image, border_color='white'):
    
    """
    detect and crop the border (white or black) of an image

    params:
    - image: NumPy array, shape should be (height, width, 3)
    - border_color: str indicating color 'white'(default) or 'black'

    return:
    - cropped image(numpy array)
    """
    
    # check if the image is RGB (3-channel)
    assert image.ndim == 3 and image.shape[2] == 3, "image must be (height, width, 3) RGB image"
    
    # set the color to crop
    if border_color in ('white', 'w'):
        color_to_crop = [255, 255, 255]
    elif border_color in ('black', 'b'):
        color_to_crop = [0, 0, 0]
    else:
        raise ValueError("only 'white' or 'black' are supported for border_color")

    # find which pixels are not the border color
    mask = np.all(image != color_to_crop, axis=-1)

    # get the border of the non-border color
    coords = np.argwhere(mask)

    # return the original image if no border is found
    if coords.shape[0] == 0:
        return image

    # get the bounding box of the non-border color
    y_min, x_min = coords.min(axis=0)
    y_max, x_max = coords.max(axis=0)

    # crop the image
    cropped_image = image[y_min:y_max+1, x_min:x_max+1]

    return cropped_image


# ----------------------------------------------------------------------
# >>> crop border 2                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/11/15  - created
#
# Desc
#
# ----------------------------------------------------------------------

def crop_to_rect_map(image, buff=0):
    """
    Crop an image to keep only the inner rectangular contour map, removing any elements outside it.

    Parameters:
    - image: NumPy array with shape (height, width, 3), an RGB image.
    - buff:  int, width of the buffer(in the number of pixels). Buffer is the region 
             where the edge of the image should falls in. 
    Returns:
    - Cropped image as a NumPy array.
    """
    
    # Convert image to grayscale for simpler processing
    gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

    # Threshold the image to get a binary image
    # the part we want to crop is white, a value very close to 255, so we set a 
    # threshold of 250 to get a binary image
    
    bw = cv2.threshold(gray, 250, 255, cv2.THRESH_BINARY)[1]

    # Detect edges using Canny edge detection
    edges = cv2.Canny(bw, 50, 150)

    # Hough transform to detect lines
    # lines are 3D matrix, [number of lines, 1, four values (x1, y1, x2, y2)]
    # minLineLength is adjusted ad hoc here to avoid fragment lines near the edges.
    
    lines = cv2.HoughLinesP(edges, 1, np.pi / 180, threshold=100, minLineLength=200, maxLineGap=10)

    # crop the image based on the range of the lines
    # image pixels order is from top-left to bottom-right (n_y, n_x)
    
    n_y, n_x = image.shape[:2]
    
    if lines is not None:
        
        xmin = min([line[0,0] for line in lines])
        ymin = min([line[0,1] for line in lines])
        xmax = max([line[0,2] for line in lines])
        ymax = max([line[0,3] for line in lines])
        
        # draw the rectangle on edges (in case you want to visualize the cropping box)
        # cv2.line(edges, (xmin, ymin), (xmin, ymax), (155, 155, 155), 2)
        # cv2.line(edges, (xmin, ymax), (xmax, ymax), (155, 155, 155), 2)
        # cv2.line(edges, (xmax, ymax), (xmax, ymin), (155, 155, 155), 2)
        # cv2.line(edges, (xmax, ymin), (xmin, ymin), (155, 155, 155), 2)

        # crop the image
        
        if buff > 0 :
            
            if n_y - ymax > buff:  ymax = n_y
            if ymin       > buff:  ymin = 0
            if n_x - xmax > buff:  xmax = n_x
            if xmin       > buff:  xmin = 0 

        cropped_image = image[ymin:ymax, xmin:xmax]

    else:
        cropped_image = image

    # return the cropped image
    return cropped_image


# ----------------------------------------------------------------------
# >>> oblique shock beta                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/03/26  - created
#
# Desc
#    - calculate the oblique shock angle beta (in degrees) given the incoming Mach 
#     M1 and flow deflection angle theta.
# ----------------------------------------------------------------------

def oblique_shock_beta(M1, theta_deg, gamma=1.4):
    
    """
    Calculate the oblique shock angle beta (in degrees) given the incoming Mach 
    number M1 and the flow deflection angle theta.
    Parameters:
        M1 (float): incoming Mach > 1
        theta_deg (float): deflection angle 0 < theta < theta_max）。
        gamma (float): specific heat ratio (default: 1.4)

    return:
        tuple: (weak_solution_deg, strong_solution_deg) weak and strong solutions
        of the oblique shock angle in degrees.
        Return None if there is no solution.
    """
    theta = np.radians(theta_deg)
    
    def theta_beta_equation(beta_rad):
        sin_beta    = np.sin(beta_rad)
        numerator   = (M1**2 * sin_beta**2 - 1)
        denominator = M1**2 * (gamma + np.cos(2 * beta_rad)) + 2.0
        return np.tan(theta) - 2 * (1 / np.tan(beta_rad)) * numerator / denominator

    # initial guess
    beta_guess_weak   = np.arcsin(1 / M1)*1.1  # weak solution guess（> close to Mach angle）
    beta_guess_strong = np.radians(80)     # strong solution guess（close to normal shock）

    # weak soluion
    sol_weak  = root(theta_beta_equation, beta_guess_weak)
    
    beta_weak = sol_weak.x[0] if sol_weak.success else None

    # strong solution
    sol_strong = root(theta_beta_equation, beta_guess_strong)
    
    beta_strong = sol_strong.x[0] if sol_strong.success else None

    # check if the solution （β > μ and θ < θ_max）
    if beta_weak is not None:
        beta_weak_deg = np.degrees(beta_weak)
        # check if the deflection angle is within the maximum deflection angle
        theta_max = oblique_shock_theta_max(M1, gamma)
        if theta_deg > theta_max:
            return None, None
    else:
        beta_weak_deg = None

    if beta_strong is not None:
        beta_strong_deg = np.degrees(beta_strong)
    else:
        beta_strong_deg = None

    return beta_weak_deg, beta_strong_deg


# ----------------------------------------------------------------------
# >>> oblique_shock_theta_max                                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/03/26  - created
#
# Desc
#  - calculate the maximum deflection angle, given the incoming Mach number M1.
# ----------------------------------------------------------------------

def oblique_shock_theta_max(M1, gamma=1.4):
    from scipy.optimize import root_scalar
    """
    Calculate the maximum flow deflection angle theta_max (in degrees) given the
    incoming Mach number M1.
    """
    def theta_beta_equation(beta_rad):
        """θ-β-Mach relations"""
        sin_beta = np.sin(beta_rad)
        numerator = M1**2 * sin_beta**2 - 1
        denominator = M1**2 * (gamma + np.cos(2 * beta_rad)) + 2
        return np.arctan(2 * numerator / (np.tan(beta_rad) * denominator))

    def dtheta_dbeta(beta_rad):
        """numerics to calculate dθ/dβ"""
        h = 1e-6
        return (theta_beta_equation(beta_rad + h) - theta_beta_equation(beta_rad - h)) / (2 * h)

    # find β（critical shock angle）where dθ/dβ = 0 
    try:
        sol = root_scalar(
            dtheta_dbeta,
            bracket=[np.arcsin(1/M1) + 0.01, np.pi/2 - 0.01],  # between Mach angle and 90°
            method='brentq'
        )
        beta_critical = sol.root
    except ValueError:
        raise ValueError("cannot find the critical shock angle, M1 too low or incorrect gamma.")
    
    # calculate the maximum deflection angle
    theta_max_rad = theta_beta_equation(beta_critical)
    return np.degrees(theta_max_rad)


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/03/26  - created
#
# Desc
#
# ----------------------------------------------------------------------

def quantities_after_oblique_shock(M1, beta, theta, gamma=1.4):
    
    """
    calculate the flow quantities ratio after an oblique shock wave.
    Parameters:
        M1 (float): incoming Mach number
        beta (float): oblique shock angle in degrees
        theta (float): flow deflection angle in degrees
        gamma (float): specific heat ratio (default: 1.4)
    Returns:
        tuple: (M2, P2/P1, T2/T1, rho2/rho1) after the oblique shock wave.
    """
    # - from J.D.Anderson, Fundamentals of Aerodynamics, 5th edition, p. 611
    
    mach_1n   = M1 * np.sin(np.radians(beta))
    mach_2n   = np.sqrt( (mach_1n**2 * (gamma-1.0)/2 + 1.0) / 
                         (gamma*mach_1n**2 - (gamma-1.0)/2) )
    rho2_rho1 = (gamma+1.0)*mach_1n**2 / (2.0 + (gamma-1.0)*mach_1n**2) 
    p2_p1     = 1.0 + 2.0*gamma/(gamma+1.0)*(mach_1n**2 - 1.0)
    t2_t1     = p2_p1 / rho2_rho1
    M2        = mach_2n / np.sin(np.radians(beta-theta))
    
    return M2, p2_p1, t2_t1, rho2_rho1 


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
    
# --- test crop_to_rect_map
    
#     imagefile = "/home/wencanwu/Downloads/figs/slice_z_01927735.png"
    
#     imag = cv2.imread(imagefile)
#     imag = cv2.cvtColor(imag, cv2.COLOR_BGR2RGB)
# #    imag = np.array(imag)
    
#     image = crop_to_rect_map(imag)
    
#     print(type( image ))
    
#     plt.imshow(image)
    
#     plt.show()

# --- test oblique_shock_beta

    M1 = 2.0
    theta_deg = 5.992
    
    
    weak, strong = oblique_shock_beta(M1, theta_deg)
    
    print(f"weak solution  : {weak:.2f} degrees")
    print(f"strong solution: {strong:.2f} degrees")
    
    print( quantities_after_oblique_shock(M1, 40.04, 0.0) )
    