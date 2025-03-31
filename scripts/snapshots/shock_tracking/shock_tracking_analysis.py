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

from   vista.params       import Params
from   vista.directories  import Directories
from   vista.psd          import pre_multi_psd
from   vista.tools        import get_filelist
from   vista.directories  import create_folder
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(fontsize=20)
# =============================================================================

casedir = '/home/wencan/temp/smooth_mid/'

fs      = 100
n_seg   = 3
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

    x_spanave = list()
    x_mid     = list()
    x_shocks  = list()
    
    for shockline in shocklines:
        x_shock = np.array( shockline['x'] )
        x_spanave.append( x_shock.mean() )
        x_mid.append( x_shock[int(len(x_shock)/2)] )
        x_shocks.append( x_shock )
        
    write_shock_loc_rms( x_spanave, x_mid )
    
    plot_shock_loc( times, x_spanave, x_mid )
    print("shock_pmpsd_spanave.png is saved.")
    
    shock_motion_psd( x_spanave )
    
    u_spanave  = np.diff( x_spanave ) / np.diff( times ) / 5.2
    u_mid      = np.diff( x_mid )     / np.diff( times ) / 5.2
    u_relative = u_mid - u_spanave
    
    u_rel_pos  = u_relative[ u_relative > 0 ]
    u_rel_neg  = u_relative[ u_relative < 0 ]
    
    times_u = times[:-1] + 0.5*np.diff(times)
    
    print(f"Mean positive relative shock velocity: {np.mean(u_rel_pos)}")
    print(f"Mean negative relative shock velocity: {np.mean(u_rel_neg)}")
    print(f"Percentage of positive relative shock velocity: {len(u_rel_pos)/len(u_relative)}")
    print(f"Percentage of negative relative shock velocity: {len(u_rel_neg)/len(u_relative)}\n")
    
    print("shock speed at the mid plane u_mid:" )
    shock_velocity_analysis( u_mid, 'u_mid_pdf.png' )
    
    print("shock speed at the spanwise average u_spanave:" )
    shock_velocity_analysis( u_spanave, 'u_spanave_pdf.png' )
    
    print("shock speed of relative u_relative (u_mid-u_span):")
    shock_velocity_analysis( u_relative, 'u_relative_pdf.png' )
    
    print("============")



def write_shock_loc_rms( x_spanave, x_mid ):
    
    x_mean = np.mean(x_spanave)
    x_fluc = x_spanave - x_mean
    rms3d  = np.sqrt( np.mean( x_fluc**2 ) )
    rms2d  = np.sqrt( np.mean( np.array(x_mid-x_mean)**2 ) )
    with open('shock_rms.dat','w') as f:
        f.write(f"mean shock location               : {x_mean}\n")
        f.write(f"spanwise averaged shock motion rms: {rms3d }\n")
        f.write(f"mid-span shock motion rms         : {rms2d }\n")
        f.write("\n=====================================\n")
        f.write("Note: values are in the original unit.\n")

def plot_shock_loc( times, x_spanave, x_mid ):

    x_mean = np.mean(x_spanave)
    
    # - plot shock motion
    
    fig, ax = plt.subplots( figsize=(12,6) )
    
    times = ( np.array(times) - 20.0 ) * 507.0 /5.2
    
    ax.plot( times, np.array(x_spanave - x_mean)/5.2,'b', label='spanwise averaged' )
    ax.plot( times, np.array(x_mid-x_mean)/5.2,'r', ls=':',label='mid-span' )
    ax.set_ylim( -0.5, 0.5 )
    ax.set_xlabel(r"$t u_{\infty}/\delta_0$")
    ax.set_ylabel(r"$x_{shock}^{'}/\delta_0$")
    ax.set_title('Shock location')
    ax.legend()
 
    plt.savefig( 'shock_location_3d.png' )
    plt.close()
    

def shock_motion_psd( x_shocks ):
    # - compute the psd of mean shock fluctuation
    
    freq, psd = pre_multi_psd( x_shocks, fs, n_seg, overlap, nfft=len(x_shocks) )
    freq = freq * Lsep / 507.0
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(freq, psd, label='3d shock tracking')
    ax.set_xscale('log')
    ax.set_xlabel(r'$St_{L_{sep}}$')
    ax.set_ylabel('pmPSD')
    
    plt.savefig( 'shock_pmpsd_spanave.png' )
    plt.close()  


def plot_spanwise_velocity_skewness( x_shocks, times ):
    
    # - plot spanwise distribution of shock velocity skewness

    x_shocks = np.array( x_shocks ).T
    
    fig, ax = plt.subplots( figsize=(12,6) )
    for i in range( len(x_shocks) ):
        print(i)
        sk = shock_velocity_analysis( times, x_shocks[i] )
        ax.plot(i*0.216667/5.2-2.0, sk, 'o', color='b')
    ax.set_ylim( -0.5, 0.0 )
    ax.set_xlabel(r"$z/\delta_0$")
    ax.set_ylabel('Skewness of shock velocity')
    plt.show()
    plt.close()


def shock_velocity_analysis( u_shock, figname ):
    
    # - compute the velocity of the shock

    skewness = skewness_analysis( u_shock, figname )
    flatness = kurtosis(u_shock)
    
    print(f"Shock velocity skewness: {skewness}, flatness: {flatness}")

    # - get all the positive shock velocity
    
    u_positive = u_shock[ u_shock > 0 ]
    u_negative = u_shock[ u_shock < 0 ]
    
    print(f"Mean positive shock velocity: {np.mean(u_positive)}")
    print(f"Mean negative shock velocity: {np.mean(u_negative)}")
    print(f"Percentage of positive shock velocity: {len(u_positive)/len(u_shock)}")
    print(f"Percentage of negative shock velocity: {len(u_negative)/len(u_shock)}")

    # fig, ax = plt.subplots( figsize=(12,6) )
    
    # ax.plot( times, u_shock, 'b', lw=0.5 )
    # ax.set_title('Shock velocity')
#    plt.show()
#    plt.close()
    
    print(f"")
    
    return skewness


def skewness_analysis( var, figname=None ):
    
    # - compute the std of the variable
    
#    var = np.array(var) - np.mean(var)
    
    std = np.std(var)
    print( f"Standard deviation of the variable: {std}" )
    print( f"range of var/std: {min(var/std)} - {max(var/std)}" )
    
    # - compute the skewness of the variable
    
    skewness = skew(var)
    print( f"Skewness of the variable: {skewness}" )
    
    # - compute the p.d.f of the variable
    
    hist, bin_edges = np.histogram( var/std, bins=80, range=(-4,4), density=True )
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:]) 
    
    # - normal distribution
    
    x = np.linspace(-5,5,100)
    y = 1/np.sqrt(2*np.pi) * np.exp(-0.5*x**2)
    
    
    fig, ax = plt.subplots( figsize=(8,6) )
    
    ax.plot( bin_centers, hist, 'o', markersize=10 )
    ax.plot( x, y, 'gray', ':', lw=1.5 )
    
    ax.set_xlabel(r"$\sigma_{u_{shock}}$")
    ax.set_ylabel('p.d.f')
    
    ax.set_xlim( -3.5, 3.5 )
    ax.set_ylim( 0, 0.5 )
    # range of the plot in the fig:
    plt.subplots_adjust(left=0.2, right=0.95, top=0.7, bottom=0.3)
    
    ax.set_yticks(np.arange(0, 0.6, 0.1))
    ax.set_yticklabels(['0.0', '0.1', '0.2', '0.3', '0.4','0.5'])
    
    ax.set_title(f'p.d.f of the variable, skewness={skewness:.2f}',pad=30)
    
    if figname is not None:
        plt.savefig( figname )
    else:
        plt.show()
    plt.close()
    
    return skewness



# =============================================================================

if __name__ == '__main__':
    
    main()