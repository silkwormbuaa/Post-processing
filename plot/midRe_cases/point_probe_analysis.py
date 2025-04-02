#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   point_probe_analysis.py
@Time    :   2025/04/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Analyse the point probe data
'''


import os
import sys
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt
from   scipy.stats       import skew, kurtosis

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.params      import Params
from   vista.directories import Directories


def main():
    casefolder = '/home/wencan/temp/smooth_mid'
    show_case(casefolder,'u')

def show_case( casefolder, var ):
    
    dirs = Directories(casefolder)
    
    point_probe_file = dirs.pp_snp_pfmax + '/point_probe.dat'
    if not os.path.exists(point_probe_file):
        print(f'File {point_probe_file} does not exist.')
        return
    
    df = pd.read_csv(point_probe_file, delimiter=r'\s+')
    
    var_fluc = df[var] - df[var].mean()
    
    rs = np.std(var_fluc)**2
    print( rs )
    
    u_fluc = df['u'] - df['u'].mean()
    v_fluc = df['v'] - df['v'].mean()
    w_fluc = df['w'] - df['w'].mean()
    p_fluc  = df['p'] - df['p'].mean()
    tke = 0.5 * (u_fluc**2 + v_fluc**2 + w_fluc**2)
    print( tke.mean() )

    
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    print( skew(var_fluc) )
    ax.plot(df['itime'], var_fluc)
    ax.set_xlabel('Time')
    ax.set_ylabel(var)
    ax.set_xlim(20,25)
    plt.show()

# =============================================================================
if __name__ == "__main__":

    main()

