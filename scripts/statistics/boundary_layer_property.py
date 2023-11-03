#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   boundary_layer_property.py
@Time    :   2023/10/28 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Compute boundary layer properties from profile 
'''


import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.line        import ProfileData

path0 = '/media/wencanwu/Seagate Expansion Drive/temp/smooth_wall/results/profile'
path1 = '/media/wencanwu/Seagate Expansion Drive/temp/221014/results/profile'
path2 = '/media/wencanwu/Seagate Expansion Drive/temp/220926/results/profile'
path3 = '/media/wencanwu/Seagate Expansion Drive/temp/220825/results/profile'
path4 = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/profile'
path5 = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results/profile'

paths = [path0, path1, path2, path3, path4, path5]
dy       = [0,  0.494, 0.468, 0.416, 0.312, 0.26]
u_ref = 507.0

# =============================================================================

for i, path in enumerate( paths ):
    
    os.chdir( path )
    
    line = ProfileData( 'profile_mean.dat' )
    line.shift_y( dy[i] )
    
    delta99,delta_star,theta = line.compute_boundary_layer_property(u_ref)
    
    h_factor = delta_star / theta
    
    with open("boundary_layer_property.dat", "w") as f:
        f.write(f"delta99    = {delta99:10.5f}\n")
        f.write(f"delta_star = {delta_star:10.5f}\n")
        f.write(f"theta      = {theta:10.5f}\n")
        f.write(f"h_factor   = {h_factor:10.5f}\n")
        
    print(f"{i} {delta99:10.5f}{delta_star:10.5f}{theta:10.5f}{h_factor:10.5f}")