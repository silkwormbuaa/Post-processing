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
from   vista.log          import Logger
from   vista.psd          import pre_multi_psd
from   vista.math_opr     import find_parabola_max
from   vista.directories  import create_folder
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(fontsize=10)
# =============================================================================

def main():
    
    case_dir = '/home/wencan/temp/241030/'
    axis     = True
    
    dirs     = Directories( case_dir )
    params   = Params( dirs.case_para_file )
    
    os.chdir( create_folder(dirs.pp_correlation))
    Logger( 'correlation_analysis.log' )
    
    dfsk1 = read_shock_loc( dirs.pp_shock + '/group1/shock_tracking1.pkl' )
    dfsk2 = read_shock_loc( dirs.pp_shock + '/group2/shock_tracking2.pkl' )
    dfbb  = read_bubble_size( dirs.pp_bubble + '/bubble_size.dat' )
    dfprb = read_prb_pressure( dirs.fetch_prb_from_type('pfmax'), params )
    dfp   = read_p_spanave_pfmax( dirs.pp_snp_pfmax + '/pressure_at_pfmax.dat' )
    
    times = ( dfsk1['itime'] - 20.0 ) * 507.0 / 5.2
    
    plot_bubble_size_fluc( times, dfbb['fluc_thr'],        params, figname='fluc_bb.png',     showaxis=axis )
    plot_bubble_size_fluc( times, dfbb['fluc'],            params, figname='fluc_bb_all.png', showaxis=axis )
    plot_shock_loc(        times, dfsk1['x_fluc_spanave'], params, figname='fluc_sk1.png', showaxis=axis )
    plot_shock_loc(        times, dfsk2['x_fluc_spanave'], params, figname='fluc_sk2.png', showaxis=axis )
    plot_pressure_fluc(    times, dfp['p_spanave_pfmax'],  params, figname='fluc_pw.png',  showaxis=axis )
    
    # print( dfbb.iloc[:10] )
    
    do_psd( dfbb['fluc_thr'],        params, figname='psd_bb.png',     showaxis=axis )
    do_psd( dfbb['fluc'],            params, figname='psd_bb_all.png', showaxis=axis )
    do_psd( dfsk1['x_fluc_spanave'], params, figname='psd_sk.png',  showaxis=axis )
    do_psd( dfp['p_fluc'],           params, figname='psd_pw.png',  showaxis=axis )
    do_psd( dfsk2['x_fluc_spanave'], params, figname='psd_sk2.png', showaxis=axis )
    
    # plot_fluc( times, dfsk1['x_fluc_spanave'], 'shock1 spanave' )
    # plot_fluc( times, dfp['p_fluc'],           'pressure at pfmax' )
    # plot_fluc( times, dfbb['fluc_thr'],        'bubble size fluctuation(y>0)' )
    # plot_fluc( times, dfbb['fluc'],            'bubble size fluctuation' )
    # plot_fluc( times, dfprb['p_fluc'],         'pressure fluctuation' )
    
    corr( dfp['p_fluc'],           dfbb['fluc_thr'], 'p_vs_bb' )
    corr( dfsk1['x_fluc_spanave'], dfp ['p_fluc'],   'shock motion vs pressure' )
#    corr( dfp['p_fluc'],           dfbb['fluc'],     'pressure vs bubble size(all)' )
    corr( dfsk1['x_fluc_spanave'], dfbb['fluc_thr'], 'shock motion vs bubble size(y>0)' )
#    corr( dfsk1['x_fluc_spanave'], dfsk2['x_fluc_spanave'],'shock motion vs shock motion2' )
#    corr( dfsk1['x_fluc_spanave'], dfbb['fluc'],     'shock motion vs bubble size(all)' )
    
#    corr( dfbb['fluc'], dfprb['p_fluc'], 'bubble size vs pressure' )
#    corr( dfsk1['x_fluc_spanave'], dfprb['p_fluc'], 'shock1 spanave vs pressure' )
    
    corr( dfsk1['x_fluc_spanave'], dfsk2['x_fluc_spanave'], 'shock at 2 delta vs shock at 6 delta' )
    

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

def read_p_spanave_pfmax( file ):
    
    # - read in the extracted spanwise averaged pressure at pfmax location
    # ['itstep', 'itime', 'p_spanave']

    df = pd.read_csv( file, delimiter=r'\s+' )

    df['p_fluc'] = df['p_spanave_pfmax'] - np.mean(df['p_spanave_pfmax'])

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

    print(f"lag at max correlation: {lag_at_max_corr}, max correlation: {max_corr} of {title}")
    
    plot_corr( lags, correlation, lag_at_max_corr, title )
    
    return lags, correlation, lag_at_max_corr


def plot_corr( lags, correlation, lag_at_max_corr, title='' ):

    fig, axs = plt.subplots( 1,2, figsize=(8,6) )
    axs[0].plot( lags, correlation )
    axs[0].axvline( lag_at_max_corr, color='r', ls='--' )
    # axs[0].set_title( title )
    
    
    axs[1].plot( lags, correlation )
    axs[1].axvline( lag_at_max_corr, color='r', ls='--' )
    axs[1].set_xlim([-50,50])
    # axs[1].set_title( title )
    
    
    plt.show()
    plt.close()
    

def plot_fluc( times, fluc, title='' ):
    
    fig, ax = plt.subplots( figsize=(8,4) )

    norm = 0.5*(np.array(fluc).max() - np.array(fluc).min())
    ax.plot( times, np.array(fluc)/norm,'b' )

    ax.set_ylim( -1.1, 1.1 )
    ax.set_title(title)

    plt.show()
    plt.close()

def adjust_fluctuation_plots( ax: plt.Axes ):

    # ax.set_title('bubble size fluctuation')

    ax.tick_params(axis='both', which='major', length=5,   width=1)
    ax.tick_params(axis='both', which='minor', length=2.5, width=1)
    ax.tick_params(axis='both', which='both',  direction='in', zorder=10, pad=5)
    ax.xaxis.set_tick_params(width=1)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines[:].set_linewidth(1)
    ax.spines[:].set_zorder(10)

def plot_bubble_size_fluc(times, fluc, params:Params, showaxis=True, figname=None):
    
    fig, ax = plt.subplots( figsize=(5,1.5) ,constrained_layout=True )

    fluc = np.array(fluc) - np.mean(fluc)
    
    norm = params.delta_0 * params.delta_0 * 20.8
    
    times = np.array( times )
    
    ax.plot( np.array([times[0],times[-1]]),  np.array([0.0,0.0]), 
            'gray', ls='--', linewidth=1.0, zorder=1 )
    ax.plot( times, np.array(fluc)/norm,'black', linewidth=1.0 )

    ax.set_xlim( 0, times[-1] )

    if showaxis:
        ax.set_xlabel(r"$(t-t_0)u_{\infty}/\delta_0$")
        ax.set_ylim( -1.1, 1.0)
    else:
        ax.spines['bottom'].set_visible(False)
        ax.set_ylim( -1.0, 1.0 )
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        
    ax.set_ylabel(r"$V_{rev}-\langle V_{rev} \rangle / \delta_0^2 L_z$")

    adjust_fluctuation_plots( ax )

    if figname is not None:
        plt.savefig( figname, dpi=600 )
    plt.show()
    plt.close()


def plot_shock_loc(times, fluc, params:Params, showaxis=True, figname=None):
    
    fig, ax = plt.subplots( figsize=(5,1.5) ,constrained_layout=True )

    fluc = np.array(fluc) - np.mean(fluc)
    
    norm = params.delta_0
    
    times = np.array( times )
    
    ax.plot( np.array([times[0],times[-1]]),  np.array([0.0,0.0]), 
            'gray', ls='--', linewidth=1.0, zorder=1 )
    ax.plot( times, np.array(fluc)/norm,'black', linewidth=1.0 )

    ax.set_xlim( 0, times[-1] )
    # ax.set_title('bubble size fluctuation')

    if showaxis:
        ax.set_xlabel(r"$(t-t_0)u_{\infty}/\delta_0$")
        ax.set_ylim( -0.55, 0.5)
    else:
        ax.spines['bottom'].set_visible(False)
        ax.set_ylim( -0.5, 0.5 )
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        
    ax.set_ylabel(r"$\langle x_{sw} - \langle x_{sw}\rangle \rangle_{z}/\delta_0$")

    adjust_fluctuation_plots( ax )

    if figname is not None:
        plt.savefig( figname, dpi=600 )
    plt.show()
    plt.close()

def plot_pressure_fluc(times, fluc, params:Params, showaxis=True, figname=None):

    fig, ax = plt.subplots( figsize=(5,1.5) ,constrained_layout=True )

    fluc = np.array(fluc) - np.mean(fluc)
    
    norm = params.p_ref
    
    times = np.array( times )
    
    ax.plot( np.array([times[0],times[-1]]),  np.array([0.0,0.0]), 
            'gray', ls='--', linewidth=1.0, zorder=1  )
    ax.plot( times, np.array(fluc)/norm,'black', linewidth=1.0 )

    ax.set_xlim( 0, times[-1] )
    # ax.set_title('bubble size fluctuation')

    if showaxis:
        ax.set_xlabel(r"$(t-t_0)u_{\infty}/\delta_0$")
        ax.set_ylim( -0.25, 0.20)
    else:
        ax.spines['bottom'].set_visible(False)
        ax.set_ylim( -0.20, 0.20 )
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    
    ax.set_ylabel(r"$\langle p_w - \langle p_w \rangle \rangle_z/p_{\infty}$")

    adjust_fluctuation_plots( ax )
    
    if figname is not None:
        plt.savefig( figname, dpi=600 )
    plt.show()
    plt.close()

def do_psd( fluc, params:Params, showaxis=True, figname=None ):
    
    freq, pm_psd = pre_multi_psd( fluc, 100, 8, 0.5, nfft=len(fluc) )
    
    freq = freq * 9.22 * params.delta_0 / params.u_ref
    
    fig, ax = plt.subplots( figsize=(2,1.5), constrained_layout=True )
    ax.semilogx( freq, pm_psd, color='black', linestyle='-', linewidth=1.0 ) 
    ax.fill_between( freq, np.zeros_like(pm_psd), pm_psd, color='gray', alpha=0.1 )
    
    if showaxis:
        ax.set_xlabel(r"$St_{L_{sep}}$")
        ax.set_ylim( [-0.05,0.9] )
    else:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        ax.set_ylim( [0.0,0.9] )
    
    ax.set_yticks( np.array([0.0, 0.5]))
    ax.set_ylabel(r"$f \cdot \mathcal{P}(f)/ \int \mathcal{P}(f) \mathrm{d} f$")
    ax.set_xlim( [1e-3, 1] )
    adjust_fluctuation_plots( ax )
    
    if figname is not None:
        plt.savefig( figname, dpi=600 )
    plt.show()
    plt.close()

# =============================================================================
if __name__ == "__main__":
    main()