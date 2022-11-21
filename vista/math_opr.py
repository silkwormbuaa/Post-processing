#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_math.py
@Time    :   2022/10/13 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Basic math functions 
'''

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
