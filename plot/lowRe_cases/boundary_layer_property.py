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

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line        import ProfileData

outpath = '/media/wencan/Expansion/temp/DataPost/lowRe/profile_arrays'

cases = ['smooth_adiabatic', '221014', '220926', '220825', '220927',  '221221']

datapaths = [f'/media/wencan/Expansion/temp/{case}/postprocess/statistics/profile_array' for case in cases]

dy       = [0,  0.494, 0.468, 0.416, 0.312, 0.26]
u_ref = 507.0

# =============================================================================

for i, path in enumerate( datapaths ):
    
    os.chdir( path )
    
    line = ProfileData( 'profile_mean_05.dat' )
    line.shift_y( dy[i] )
    
    delta99,delta_star,theta = line.compute_boundary_layer_property(u_ref)
    
    h_factor = delta_star / theta
    
    with open("boundary_layer_property_05.dat", "w") as f:
        f.write(f"delta99    = {delta99:10.5f}\n")
        f.write(f"delta_star = {delta_star:10.5f}\n")
        f.write(f"theta      = {theta:10.5f}\n")
        f.write(f"h_factor   = {h_factor:10.5f}\n")
        
    print(f"{i} {delta99:10.5f}{delta_star:10.5f}{theta:10.5f}{h_factor:10.5f}")