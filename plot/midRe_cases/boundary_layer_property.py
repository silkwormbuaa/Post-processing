#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   boundary_layer_property.py
@Time    :   2025/02/03 
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

cases     = ['smooth_mid', '231124', '241030','smooth_adiabatic','220927'] #
datapaths = [f'/media/wencan/Expansion/temp/{case}/postprocess/statistics/profile_incip' for case in cases]
dy        = [0, 0.312, 0.07921,0.0,0.312]
u_ref     = 507.0

# =============================================================================

for i, path in enumerate( datapaths ):
    
    os.chdir( path )
    
    line = ProfileData( 'profile_mean.dat' )
    line.shift_y( dy[i] )
    
    # drop some values at the tail
    line.df = line.df.iloc[:-5]
    line.df.reset_index(drop=True, inplace=True)
    
    delta99,delta_star,theta = line.compute_boundary_layer_property(u_ref)
    
    h_factor = delta_star / theta
    
    with open("boundary_layer_property_incip.dat", "w") as f:
        f.write(f"delta99    = {delta99:10.5f}\n")
        f.write(f"delta_star = {delta_star:10.5f}\n")
        f.write(f"theta      = {theta:10.5f}\n")
        f.write(f"h_factor   = {h_factor:10.5f}\n")
        
    print(f"{i} {delta99:10.5f}{delta_star:10.5f}{theta:10.5f}{h_factor:10.5f}")