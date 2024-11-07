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

shockpath2d = casedir + '/postprocess/snapshots/shock_tracking/2d_y2'
shockpath3d = casedir + '/postprocess/snapshots/shock_tracking/3d'

fs = 100
n_seg = 8
overlap = 0.5

analyse2d=True
analyse3d=False

dirs = Directories( casedir )

# =============================================================================

params = Params( dirs.case_para_file )
Lsep   = params.Lsep


if analyse2d:
    shock2d_files = get_filelist( shockpath2d, 'shock_tracking2d' )
if analyse3d:
    shock3d_files = get_filelist( shockpath3d, 'shock_tracking3d' )

if analyse2d:
    
    os.chdir( create_folder( shockpath2d ) )
    
    # - read the shock motion data into dataframe    
    
    for shock2d_file in shock2d_files:
        df2d = pd.concat([pd.read_csv(shock2d_file, delimiter=r"\s+") for shock2d_file in shock2d_files])
    
    # - compute the rms of the shock fluctuation
    
    x_mean = df2d['x_shock'].mean()
    x_fluc = df2d['x_shock'] - x_mean
    rms2d  = np.sqrt( np.mean( x_fluc**2 ) )
    skewness = skew( x_fluc )
    flatness = kurtosis( x_fluc )
    print(f" rms: {rms2d}, skewness: {skewness}, flatness: {flatness}")
    
    # - plot shock motion
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot( df2d['time'], df2d['x_shock'] )
    ax.set_title('Shock location at Z=0.0 plane')

    plt.savefig( 'shock_location_2d.png' )
    plt.close()
        
    # - compute the psd of the shock fluctuation
    
    f2d, psd2d = pre_multi_psd( x_fluc, fs, n_seg, overlap, nfft=len(x_fluc) )
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(f2d, psd2d, label='2d shock tracking')
    
    ax.set_xscale('log')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('pmPSD')

    plt.savefig( 'shock_psd_2d.png' )
    plt.close()

if analyse3d:
    
    os.chdir( create_folder( shockpath3d ) )
    
    # - read in the shock motion data into dataframe
    
    times = list()
    shocklines = list()
    
    for shock3d_file in shock3d_files:
        with open(shock3d_file, 'rb') as f:
            time = pickle.load(f) 
            shockline = pickle.load(f)

        times = times + time
        shocklines = shocklines + shockline
    
    print(len(times),len(shocklines))
    
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
    