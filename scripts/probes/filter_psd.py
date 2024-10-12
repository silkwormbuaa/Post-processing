#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   filter_psd.py
@Time    :   2024/10/12 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Compute psd of probe signals, filter and reconstruct the signals
'''

import os
import sys
import pickle
import numpy              as     np
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.timer        import timer
from   vista.probe        import ProbeData
from   vista.probe        import ProbeFile
from   vista.directories  import Directories
from   vista.tools        import find_indices
from   vista.tools        import get_filelist
from   vista.tools        import read_case_parameter
from   vista.directories  import create_folder
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams()

# =============================================================================
# option zone
# =============================================================================

case_dir = '/home/wencanwu/test/231124/'

# =============================================================================

dirs     = Directories( case_dir )
probes   = ProbeFile( dirs.set_prb )

# -- enter working directory and read relevant parameters

os.chdir( dirs.prb_dir )

prb_files = dirs.probes
n_data = len( prb_files )
print(f"We are in directory:{dirs.prb_dir}\n")
print(f"We have got {n_data:5d} probes data.\n")

parameters = read_case_parameter( dirs.case_para_file )
p_ref   = float( parameters.get('p_ref') )
u_ref   = float( parameters.get('u_ref') )
h       = float( parameters.get('H') )
x_imp   = float( parameters.get('x_imp') )
delta_0 = float( parameters.get('delta_0') )
# the following parameters are normalized by delta_0
x_sep   = float( parameters.get('x_sep') )
x_att   = float( parameters.get('x_att') )
x_pfmax = float( parameters.get('x_pfmax') )

prb_withT = True if parameters.get('prb_withT').lower() == 'true' else False

# -- create output directories

create_folder( dirs.pp_pre_ridge )
os.chdir( dirs.pp_pre_ridge )

# -- start read pressure at the ridge

x_locs  = list()
pf_rmss = list()

if not os.path.exists( 'streamwise_rms.pkl' ):
    
    with timer("Reading pressure data"):
        
# ----- read pressure data at the ridge probe by probe

        for i in range( 337 ):
            
            probe = probes.probes[i]
            xyz = probe.xyz
            
            # find the probes at the ridges
            if not(abs(xyz[1]) < 0.001 and abs(xyz[2]) < 0.001):
                continue
            
            prb_data = ProbeData( prb_files[i], withT=prb_withT )
            prb_data.cleandata( t_start=20.0 )
            prb_data.get_fluc( 'p' )
            pf = np.array( prb_data.df['p_fluc'] )
            
            x_locs.append( (xyz[0]-x_imp)/delta_0 )
            pf_rmss.append( np.sqrt( np.mean( pf**2 ) )/p_ref )
            
            print(f"Processed pressure data at probe {i} at x={xyz[0]:.2f}.")

        # -- save data into a pickle file
        with open( 'streamwise_rms.pkl', 'wb' ) as f:
            pickle.dump( [x_locs, pf_rmss], f )
            print(f"Data saved into {dirs.pp_pre_ridge}/streamwise_rms.pkl.")
    
else:
    
    with open( 'streamwise_rms.pkl', 'rb' ) as f:
        x_locs, pf_rmss = pickle.load( f )
        print(f"Data loaded from {dirs.pp_pre_ridge}/streamwise_rms.pkl.")
        
# =============================================================================
# plot streamwise rms of pressure fluctuation
# =============================================================================

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
ax.plot( x_locs, pf_rmss, '-', color='blue', linewidth=2 )

def plot_style():
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))

    ax.minorticks_on()

    ax.tick_params( which='major',
                    axis='both',
                    direction='in',
                    length=15,
                    width=2)

    ax.tick_params( which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1)

    ax.set_xlim([-15,10])
    ax.set_ylim([0.00,0.09])
    ax.set_xlabel(r"$(x-x_{imp})/\delta_0$", labelpad=-5 )  
    ax.tick_params(axis='x', pad=15)

    ax.tick_params(axis='y', pad=10)
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)

plot_style()
ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_{\infty}$" )
figname = "pressure_fluctuation_rms"

plt.savefig( figname + '.png' )
plt.close()

# =============================================================================
# filter and reconstruct the signals
# =============================================================================

psdfilelist = get_filelist( dirs.pp_psd_ridge )

x_locs = []; rms_psds = []; powers = []; powers_lp = []; powers_hp = []

for i, psdfile in enumerate( psdfilelist ):
    
    print(f"Processing {i+1:005d} {psdfile}...")
    
    probe = ProbeData()
    probe.read_psd( psdfile )
    freq   = np.array( probe.psd_df['freq'] )
    st     = freq*(x_att-x_sep)*delta_0/u_ref
    _,i_cr = find_indices( st, 1.0 )
    
    psd_p_fluc = np.array( probe.psd_df['psd_p_fluc'] )
    power      = np.sum( psd_p_fluc )*freq[1]
    powers   .append( power )
    rms_psds .append( np.sqrt(power) )
    powers_lp.append( np.sum( psd_p_fluc[:i_cr] )*freq[1] )
    powers_hp.append( np.sum( psd_p_fluc[i_cr:] )*freq[1] )
    x_locs.append( (probe.xyz[0]-x_imp)/delta_0 )

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
ax.plot( x_locs, np.array(rms_psds)/p_ref, '-', color='blue', linewidth=2 )
plot_style()
ax.set_ylabel(r"$\sqrt{\int \mathcal{P}(f) \mathrm{d} f}/p_{\infty}$" )
figname = "pressure_fluctuation_rms_psd"
plt.savefig( figname + '.png' )
plt.close()

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
ax.plot( x_locs, np.array(powers)/p_ref**2, '-', color='blue', linewidth=2 )
plot_style()
ax.set_ylabel(r"$\int \mathcal{P}(f) \mathrm{d} f/p^2_{\infty}$" )
ax.set_ylim([0,0.006])
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.002))
figname = "pressure_fluctuation_power"
plt.savefig( figname + '.png' )
plt.close()

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
ax.plot( x_locs, np.sqrt(powers_lp)/p_ref, '-', color='blue', linewidth=2 )
plot_style()
ax.set_ylabel(r"$\sqrt{\int \mathcal{P}(f) \mathrm{d} f}_{|St<1}/p_{\infty}$" )
figname = "pressure_fluctuation_power_lp"
plt.savefig( figname + '.png' )
plt.close()

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
ax.plot( x_locs, np.sqrt(powers_hp)/p_ref, '-', color='blue', linewidth=2 )
plot_style()
ax.set_ylabel(r"$\sqrt{\int \mathcal{P}(f) \mathrm{d} f}_{|St>1}/p_{\infty}$" )
figname = "pressure_fluctuation_power_hp"
plt.savefig( figname + '.png' )
plt.close()

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
ax.plot( x_locs, np.array(powers_lp)/np.array(powers), '-', color='blue', linewidth=2 )
plot_style()
ax.set_ylabel(r"$\int \mathcal{P}(f) \mathrm{d} f_{|St<1}/\int \mathcal{P}(f) \mathrm{d} f$" )
ax.set_ylim([0,1])
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
figname = "pressure_fluctuation_power_lp_ratio"
plt.savefig( figname + '.png' )
plt.close()
