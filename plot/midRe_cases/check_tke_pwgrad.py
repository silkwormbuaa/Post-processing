#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   check_tke_pwgrad.py
@Time    :   2025/04/04 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import numpy              as     np
import pandas             as     pd
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker
import matplotlib.markers as     markers

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.line         import ProfileData
from   vista.params       import Params
from   vista.directories  import Directories
from   vista.tools        import find_indices

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 20


def main():
    
    case_dirs  =  [ 'smooth_adiabatic','221014','220926','220825','220927','221221',
                    '240210',          '240211']
    case_dirs  = ['smooth_mid',      '231124','241030','241018']
    
    tkes     = list()
    pw_grads = list()
    
    fig, ax = plt.subplots( figsize=(12, 8) )
    for case in case_dirs:
        
        tke, pw_grad = get_tke_pwgrad( '/home/wencan/temp/' + case )
        
        tkes.append( tke )
        pw_grads.append( pw_grad )
        print(f"{case} {tke:10.5f} {pw_grad:10.5f}")
    
        ax.plot( tke, pw_grad, 'o', label=case,
                 markersize=15, markerfacecolor='none')
    
    ax.legend()
    ax.set_xlabel( r"$tke_{max} \cdot 100 $" )
    ax.set_ylabel( r"$\nabla p$")
    plt.show()

        

def get_tke_pwgrad( case_dir ):
    
    dirs   = Directories( case_dir )
    params = Params( dirs.case_para_file )
    
    line = ProfileData( dirs.pp_profile_incip + '/incip_profile_mean.dat' )
    line.shift_y( params.H - params.H_md )
    line.inner_normalize( dirs.pp_profile_incip + '/incip_wall_statistics.dat' )
    line.vd_transform()
    
    line.df.to_string( dirs.pp_profile_incip + '/profile_normalized.dat',
                       index=False,
                       float_format='%15.7f',
                       justify='left' )
    
    # - simply the maximum value of the spanwise averaged tke profile
    # tke = np.array(line.df['tke']) / params.u_ref**2 * 100
    
    # - integral of the tke profile
    y = np.array( line.df['y'] )
    _, index = find_indices( y, 10.4 )
    
    tke = np.array( line.df['u`u`'] )[:index] * np.array( line.df['hy'] )[:index]
    
    tke = np.max( line.df['u`u`'] )
    
    return tke, params.pw_grad_max


def check_massflowtke_pwgrad():
    
    file  = '/home/wencan/temp/DataPost/midRe/tke_incip/tke_incip.dat'

    df    = pd.read_csv( file, delimiter=r'\s+' )
    
    fig, ax = plt.subplots( figsize=(12, 8) )
    
    for i, row in df.iterrows():
        case = row['case']
        ave  = row['ave']
        tke  = row['tke']
        pwgrad = get_pwgrad(case)
        ax.plot( ave, pwgrad, 'o', label=case,
                 markersize=15, markerfacecolor='none')
    
    ax.set_xlabel( r"tke_ave" )
    ax.set_ylabel( r"$\nabla p$")
    ax.legend( fontsize=15 )
    plt.show()

def get_pwgrad(casecode):
    
    params = Params( source_dir + 'database/parameters/' + casecode )
    
    return params.pw_grad_max

# =============================================================================
if __name__ == "__main__":

    main()
    #check_massflowtke_pwgrad()
