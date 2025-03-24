#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   shock_tracking_analysis.py
@Time    :   2024/08/23 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Analysis of the shock tracking results
'''


import os
import sys
import pickle
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt
from   scipy.stats       import skew, kurtosis

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.params      import Params
from   vista.directories import Directories
from   vista.psd         import pre_multi_psd
from   vista.tools       import get_filelist
from   vista.directories import create_folder


# =============================================================================

casedir = '/home/wencan/temp/smooth_mid'

fs      = 100
n_seg   = 8
overlap = 0.5

dirs    = Directories( casedir )
params  = Params( dirs.case_para_file )
Lsep    = params.Lsep


def main():
    
    shockpath1 = dirs.pp_shock + '/group1'
    analyze_3d_shock( shockpath1 )
    
    shockpath2 = dirs.pp_shock + '/group2'
    analyze_3d_shock( shockpath2 )
    
    print("Analysis done.")

# =============================================================================

def analyze_3d_shock( shockpath ):
    
    shock3d_file = get_filelist( shockpath, 'shock_tracking' )[0]
    
    os.chdir( create_folder( shockpath ) )
    
    # - read in the shock motion data into dataframe
    
    print(f"reading {shock3d_file}")
    
    with open(shock3d_file, 'rb') as f:
        times = pickle.load(f) 
        shocklines = pickle.load(f)

    
    x_shocks = list()
    x_shocks_mid = list()
    
    for shockline in shocklines:
        x_shock = np.array( shockline['x'] )
        x_shocks.append( x_shock.mean() )
        x_shocks_mid.append( x_shock[int(len(x_shock)/2)] )
        
    x_mean = np.mean(x_shocks)
    
    print(f"mean shock location is {x_mean}")
    
    x_fluc = x_shocks - x_mean
    rms3d  = np.sqrt( np.mean( x_fluc**2 ) )
    rms2d  = np.sqrt( np.mean( np.array(x_shocks_mid-x_mean)**2 ) )
    with open('shock_rms.dat','w') as f:
        f.write(f"spanwise averaged shock motion rms: {rms3d}\n")
        f.write(f"mid-span shock motion rms         : {rms2d}\n")
        f.write("\n=====================================\n")
        f.write("Note: values are in the original unit.\n")
    
    # - plot shock motion
    
    fig, ax = plt.subplots( figsize=(12,6) )
    
    times = ( np.array(times) - 20.0 ) * 507.0 * (1.0/5.2)
    
    ax.plot( times, np.array(x_fluc)/5.2,'b' )
    ax.plot( times, np.array(x_shocks_mid-x_mean)/5.2,'r', ls=':' )
    ax.set_ylim( -0.5, 0.5 )
    ax.set_title('Spanwise averaged shock location')
 
    plt.savefig( 'shock_location_3d.png' )
    plt.close()
    
    # - compute the psd of mean shock fluctuation
    
    f3d, psd3d = pre_multi_psd( x_shocks, fs, n_seg, overlap, nfft=len(x_shocks) )
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    
    f3d = f3d * Lsep / 507.0
    ax.plot(f3d, psd3d, label='3d shock tracking')
    
    ax.set_xscale('log')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('pmPSD')
    
    plt.savefig( 'shock_pmpsd_3d.png' )
    plt.close()  


# =============================================================================

if __name__ == '__main__':
    
    main()