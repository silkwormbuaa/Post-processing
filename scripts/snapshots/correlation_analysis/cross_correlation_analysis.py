#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   cross_corelation_analysis.py
@Time    :   2024/09/09 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Cross correlation analysis of the shock motion and bubble size
'''

import os
import sys
import pickle
import numpy              as     np
import pandas             as     pd
import matplotlib.pyplot  as     plt
from   scipy.signal       import correlate

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.params       import Params
from   vista.probe        import ProbeData
from   vista.directories  import Directories
from   vista.math_opr     import find_parabola_max
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(fontsize=30)
# =============================================================================

def main():
    
    case_dir = '/home/wencan/temp/smooth_mid/'
    dirs     = Directories( case_dir )
    params   = Params( dirs.case_para_file )
    
    dfsk1 = read_shock_loc( dirs.pp_shock + '/group1/shock_tracking1.pkl' )
    dfsk2 = read_shock_loc( dirs.pp_shock + '/group2/shock_tracking2.pkl' )
    dfbb  = read_bubble_size( dirs.pp_bubble + '/bubble_size.dat' )
    dfprb = read_prb_pressure( dirs.fetch_prb_from_type('pfmax'), params )
    
    times = ( dfsk1['itime'] - 20.0 ) * 507.0 / 5.2
    
    # plot_fluc( times, dfsk1['x_fluc_spanave'], 'shock1 spanave' )
    # plot_fluc( times, dfbb['fluc_thr'],        'bubble size fluctuation(y>0)' )
    # plot_fluc( times, dfbb['fluc'],            'bubble size fluctuation' )
    # plot_fluc( times, dfprb['p_fluc'],         'pressure fluctuation' )
    
    corr( dfsk1['x_fluc_mid'], dfbb['fluc'], 'shock1 spanave vs bubble size' )
    corr( dfsk1['x_fluc_mid'], dfbb['fluc_thr'], 'shock1 spanave vs bubble size(y>0)' )
    corr( dfsk1['x_fluc_mid'], dfprb['p_fluc'], 'shock1 spanave vs pressure' )
    corr( dfbb['fluc'], dfprb['p_fluc'], 'bubble size vs pressure' )
    
    # plot_fluc( times, dfsk2['x_fluc_spanave'], 'shock2 spanave' )
    # corr( dfsk1['x_fluc_spanave'], dfsk2['x_fluc_spanave'], 'shock at 2 delta vs shock at 6 delta' )
    

def read_shock_loc( file ):
    
    with open(file, 'rb') as f:
        times = pickle.load(f) 
        shocklines = pickle.load(f)
    
    x_spanave = list()
    x_mid     = list()

    for shockline in shocklines:
        x_shock = np.array(shockline['x'])
        x_spanave.append( np.mean(x_shock) )
        x_mid.append( x_shock[int(len(x_shock)/2)] )
    
    x_mean = np.mean(x_spanave)
    
    x_fluc_spanave = np.array(x_spanave) - x_mean
    x_fluc_mid     = np.array(x_mid)     - x_mean
    
    df = pd.DataFrame( {'itime':times, 'x_fluc_spanave':x_fluc_spanave, 'x_fluc_mid':x_fluc_mid} )
    
    return df

def read_bubble_size( file ):

    # - read in the bubble size data ['itstep', 'itime', 'bubble_volume', 'bubble_volume_thr']

    df = pd.read_csv( file, delimiter=r'\s+' )

    df['fluc']     = df['bubble_volume']     - np.mean(df['bubble_volume'])
    df['fluc_thr'] = df['bubble_volume_thr'] - np.mean(df['bubble_volume_thr'])

    return df

def read_prb_pressure( file, params: Params ):
    
    # - read in the probe data 
    
    prb_data = ProbeData( file, params.prb_withT )
    prb_data.cleandata( t_start=20.0 )
    prb_data.get_fluc( ['p'] )
    index = prb_data.time_index( np.linspace(20.0, 61.0, 4101) )
    df    = prb_data.df.iloc[index]

    return df



def corr( data1, data2, title, mode='full', method='auto' ):

    # - plot shock motion
    
    len_data = len(data1)
    
    correlation = correlate( data1, data2, mode=mode, method=method )

    correlation /= (np.std(data1) * np.std(data2) * len_data)

    lags = np.arange(-len_data+1,len_data)*0.01*507/5.2   # 0.01 is the time interval

    index = np.argmax(abs(correlation))

    p1 = [lags[index-1], correlation[index-1]]
    p2 = [lags[index  ], correlation[index  ]]
    p3 = [lags[index+1], correlation[index+1]]
    
    lag_at_max_corr, max_corr = find_parabola_max( p1, p2, p3 )

    print(f"lag at max correlation: {lag_at_max_corr}")
    
    plot_corr( lags, correlation, lag_at_max_corr, title )
    
    return lags, correlation, lag_at_max_corr


def plot_corr( lags, correlation, lag_at_max_corr, title='' ):

    fig, ax = plt.subplots( figsize=(12,6) )
    ax.plot( lags, correlation )
    ax.axvline( lag_at_max_corr, color='r', ls='--' )
    ax.set_title( title )
    ax.set_xlim([-100,100])
    plt.show()
    plt.close()
    

def plot_fluc( times, fluc, title='' ):
    
    fig, ax = plt.subplots( figsize=(15,8) )

    norm = 0.5*(np.array(fluc).max() - np.array(fluc).min())
    ax.plot( times, np.array(fluc)/norm,'b' )

    ax.set_ylim( -1.1, 1.1 )
    ax.set_title(title)

    plt.show()
    plt.close()


# =============================================================================
if __name__ == "__main__":
    main()