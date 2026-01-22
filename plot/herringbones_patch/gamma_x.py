#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   gamma_x.py
@Time    :   2026/01/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import numpy             as     np
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 24

def main():
    case = '250710'
    output_folder = f'/home/wencan/temp/{case}/postprocess/statistics/secondary_flow_intensity/'
    os.chdir( output_folder )
    
    with open( 'secondary_flow_intensity.txt', 'r' ) as f_in:
        lines   = f_in.readlines()[1:]  # skip header line
        x       = []
        x_delta = []
        gamma   = []
        gamma_w = []
        for line in lines:
            parts = line.split()
            x.append(       float(parts[0]) )
            x_delta.append( float(parts[1]) )
            gamma.append(   float(parts[2]) )
            gamma_w.append( float(parts[3]) )
    
    

    fig, ax = plt.subplots(figsize=(8,6),constrained_layout=True)
    ax.plot( x_delta[::2], np.array(gamma)  [::2]*100, 'o-', color='black', linewidth=3, markersize=10 )
    ax.plot( x_delta[::2], np.array(gamma_w)[::2]*100, 's--', color='blue', linewidth=3, markersize=10 )
    ax.set_xlabel( r'$(x-x_{imp})/\delta_0$' )
    ax.set_ylabel( r'$I\cdot 100$')
    ax.set_xlim( -23, -14 )
    ax.set_ylim( 0.0, 3.0 )
    ax.tick_params(which='major',
                   axis='both',
                   direction='in',
                   length=20,
                   width=2.0)
    ax.tick_params(which='minor',
                   axis='both', 
                   direction='in',
                   length=10,
                   width=1.5)
    
    ax.spines[:].set_linewidth(2.0)
    
    plt.savefig( 'CD2_vortex_intensity.pdf' )
    plt.show()
    plt.close()

# =============================================================================
if __name__ == "__main__":

    main()

