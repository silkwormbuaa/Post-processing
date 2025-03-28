# -*- coding: utf-8 -*-
'''
@File    :   vista_math.py
@Time    :   2022/10/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Basic math functions 
'''

import numpy             as     np

# ----------------------------------------------------------------------
# >>> LINEAR INTERPOLATION                                        ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/13  - created
#
# Input
# 
# - x ; target location
# - x1, y1 ; first data point
# - x2, y2 ; second data point
# ----------------------------------------------------------------------

def lin_interpo( x, x1, y1, x2, y2 ):
    
    slope = (y2-y1) / (x2-x1)
    
    y = slope * (x-x1) + y1
    
    return y


# ----------------------------------------------------------------------
# >>> VELOCITY GRADIENT IN Y DIRECTION                            ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/13  - created
#
# Input
#
# - u1, y1; velocity and y on the first layer of grid
# - u2, y2; '''                   second         '''
# ----------------------------------------------------------------------

def du_dy( u1, y1, u2, y2, opt=1 ):
    
    if opt == 1:            #  linear gradient
        
        dudy = abs( u1 / y1 )
    
    if opt == 2:            #  second order accuracy
        
        temp1 = u1*y2*y2 - u2*y1*y1
        
        temp2 = y1 * y2 * (y2-y1)
        
        dudy = abs( temp1 / temp2 )
        
    return dudy 


# ----------------------------------------------------------------------
# >>> Unitize                                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/16  - created
#
# Desc
#   - unitize a vector to have L2 = 1
# ----------------------------------------------------------------------

def unitize_L2( vector ):
    
    l2_norm = np.linalg.norm( vector )
    
    vector_unitized = vector / l2_norm
    
    return vector_unitized


# ----------------------------------------------------------------------
# >>> find the maximum/minimum point of a parabola curve        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/11/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

def find_parabola_max( point1, point2, point3):
    
    # get the x, y coordinates of the three points
    x1, y1 = point1
    x2, y2 = point2
    x3, y3 = point3
    
    # solve the coefficients: a, b, c
    A = np.array([
        [x1**2, x1, 1],
        [x2**2, x2, 1],
        [x3**2, x3, 1]
    ])
    B = np.array([y1, y2, y3])
    
    a, b, c = np.linalg.solve(A, B)
    
    # find the x coordinate of the maximum point
    x_max = -b / (2 * a)
    
    # get the maximum value
    y_max = a * x_max**2 + b * x_max + c
    
    return (x_max, y_max)
