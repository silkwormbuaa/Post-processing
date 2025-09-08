#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   compare_profile.py
@Time    :   2025/08/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys

import matplotlib
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker
import matplotlib.markers as     markers
import numpy              as     np
import pandas             as     pd

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line         import ProfileData

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 40

def main():
    
    locs = ['_-13','_-12.5','_-12']

    data0 = f'/media/wencan/Expansion/temp/smooth_adiabatic/postprocess/statistics/'
    data1 = f'/media/wencan/Expansion/temp/250710/postprocess/statistics/'

    lines0 = []; lines1 = []

    for i, loc in enumerate(locs):
        
        os.chdir( os.path.join(data1, loc) )
        
        line0 = ProfileData(os.path.join(data0, loc, 'profile_mean.dat'))
        line0.shift_y( 0.0 )
        line0.inner_normalize(os.path.join(data0, loc, 'wall_statistics.dat'))
        line0.vd_transform()
        
        line1 = ProfileData(os.path.join(data1, loc, 'profile_mean.dat'))
        line1.shift_y( 0.0 )
        line1.inner_normalize(os.path.join(data1, loc, 'wall_statistics.dat'))
        line1.vd_transform()
    
        lines0.append(line0)
        lines1.append(line1)
    
        fig, ax = plt.subplots( figsize=(12, 6) )
        ax.plot( lines0[-1].df['u'], lines0[0].df['y'], 'black' )
        ax.plot( lines1[-1].df['u'], lines1[0].df['y'], 'orangered' )
        
        plt.savefig('u_profile.png')
        plt.show()
        plt.close(fig)

# =============================================================================
if __name__ == "__main__":

    main()

