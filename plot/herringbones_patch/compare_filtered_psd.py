#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   compare_filtered_psd.py
@Time    :   2025/08/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import pickle
import numpy              as     np
import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker

source_dir = os.path.realpath(__file__).split('plot')[0] 
sys.path.append( source_dir )

from   vista.timer        import timer
from   vista.probe        import ProbeData
from   vista.probe        import ProbeFile
from   vista.params       import Params
from   vista.directories  import Directories
from   vista.tools        import find_indices
from   vista.tools        import get_filelist
from   vista.directories  import create_folder
from   vista.plot_setting import set_plt_rcparams
from   vista.tools        import create_linear_interpolator

set_plt_rcparams()

# =============================================================================
# option zone
# =============================================================================

cutoff_st       = 0.4
independent_len = False
addmarker       = True
cases     = ['smooth_adiabatic','250821','250821']
colors    = ['gray',  'red', 'blue']
lstyle    = ['-',     '-.',  '--'] 
locations = ['ridge', 'ridge', 'valley']
outdir    = '/home/wencan/temp/DataPost/herringbones_patch/psd'

# =============================================================================

os.chdir( outdir )

dirs = [ Directories( os.path.join( '/home/wencan/temp', case ) ) for case in cases ]

params  = [ Params( dir.case_para_file ) for dir in dirs ]

if independent_len:
    lseps   = [ param.lsep  for param in params ]
else:
    lseps   = [ 9.52 ]*len(cases)

x_sep   = [ param.x_sep   for param in params ]
x_att   = [ param.x_att   for param in params ]
x_pfmax = [ param.x_pfmax for param in params ]

x_imp   = params[0].x_imp
delta_0 = params[0].delta_0
u_ref   = params[0].u_ref
p_ref   = params[0].p_ref

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

    ax.set_xlim([-20,10])
    ax.set_ylim([0.00,0.09])
    ax.set_xlabel(r"$(x-x_{imp})/\delta_0$", labelpad=-5 )  
    ax.tick_params(axis='x', pad=15)

    ax.tick_params(axis='y', pad=10)
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)


def add_sep_pfmax(i,x,y):
    interpolator = create_linear_interpolator( x, y )
    ax.plot( x_sep[i],   interpolator(x_sep[i]),   'p', color=colors[i], ms=20)
    ax.plot( x_att[i],   interpolator(x_att[i]),   'p', color=colors[i], ms=20)
    ax.plot( x_pfmax[i], interpolator(x_pfmax[i]), '*', color=colors[i], ms=20)


# =============================================================================
# filter and reconstruct the signals
# =============================================================================

x_locs_lines = []
rms_psds_lines = []
powers_lines = []
powers_lp_lines = []
powers_hp_lines = []

for i in range( len(cases) ):

    psdfilelist = get_filelist( dirs[i].pp_probes+'/psd_'+locations[i]  )

    x_locs = []; rms_psds = []; powers = []; powers_lp = []; powers_hp = []

    for j, psdfile in enumerate( psdfilelist ):
        
        print(f"Processing {i+1:005d} {psdfile}...")
        
        probe = ProbeData()
        probe.read_psd( psdfile )
        freq   = np.array( probe.psd_df['freq'] )
        st     = freq*(lseps[i])*delta_0/u_ref
        _,i_cr = find_indices( st, cutoff_st ) 
#        _,i_cr = find_indices( freq, 10.575 )     # critical frequency St_smooth=1
        
        psd_p_fluc = np.array( probe.psd_df['psd_p_fluc'] )
        power      = np.sum( psd_p_fluc )*freq[1]
        powers   .append( power )
        rms_psds .append( np.sqrt(power) )
        powers_lp.append( np.sum( psd_p_fluc[:i_cr] )*freq[1] )
        powers_hp.append( np.sum( psd_p_fluc[i_cr:] )*freq[1] )
        x_locs.append( (probe.xyz[0]-x_imp)/delta_0 )

    x_locs_lines.append( x_locs )
    rms_psds_lines.append( rms_psds )
    powers_lines.append( powers )
    powers_lp_lines.append( powers_lp )
    powers_hp_lines.append( powers_hp )

# =============================================================================
# plot the pre-multiplied psd of pressure fluctuation

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
for i in range( len(cases) ):
    ax.plot( x_locs_lines[i], np.array(rms_psds_lines[i])/p_ref, 
             ls=lstyle[i], linewidth=4, color=colors[i] )
    if addmarker:
        add_sep_pfmax(i, x_locs_lines[i], np.array(rms_psds_lines[i])/p_ref)
    
plot_style()
ax.set_ylabel(r"$\sqrt{\int \mathcal{P}(f) \mathrm{d} f}/p_{\infty}$" )
figname = "pf_rms_psd"
plt.savefig( figname + '.png' )
plt.close()

# =============================================================================
# plot the power spectra

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
for i in range( len(cases) ):
    ax.plot( x_locs_lines[i], np.array(powers_lines[i])/p_ref**2, 
             ls=lstyle[i], linewidth=4, color=colors[i] )
    if addmarker:
        add_sep_pfmax(i, x_locs_lines[i], np.array(powers_lines[i])/p_ref**2)

plot_style()
ax.set_ylabel(r"$\int \mathcal{P}(f) \mathrm{d} f/p^2_{\infty}$" )
ax.set_ylim([0,0.008])
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.002))
figname = "pf_power"
plt.savefig( figname + '.png' )
plt.close()

# =============================================================================
# plot the low-frequency-pass power spectra

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
for i in range( len(cases) ):
    ax.plot( x_locs_lines[i], np.sqrt(powers_lp_lines[i])/p_ref, 
             ls=lstyle[i], linewidth=4, color=colors[i] )
    if addmarker:
        add_sep_pfmax(i, x_locs_lines[i], np.sqrt(powers_lp_lines[i])/p_ref)
    
plot_style()
ax.set_ylabel(r"$\sqrt{\int \mathcal{P}(f) \mathrm{d} f}_{|St<"+f"{cutoff_st}"+r"}/p_{\infty}$" )
figname = "pf_rms_lp"
plt.savefig( figname + '.png' )
plt.close()

# =============================================================================
# plot the high-frequency-pass power spectra

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
for i in range( len(cases) ):
    ax.plot( x_locs_lines[i], np.sqrt(powers_hp_lines[i])/p_ref, 
             ls=lstyle[i], linewidth=4, color=colors[i] )
    if addmarker:
        add_sep_pfmax(i, x_locs_lines[i], np.sqrt(powers_hp_lines[i])/p_ref)
    
plot_style()
ax.set_ylabel(r"$\sqrt{\int \mathcal{P}(f) \mathrm{d} f}_{|St>"+f"{cutoff_st}"+r"}/p_{\infty}$" )
figname = "pf_rms_hp"
plt.savefig( figname + '.png' )
plt.close()

# =============================================================================
# plot the ratio of low-frequency-pass power spectra to the total power spectra

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
for i in range( len(cases) ):
    ax.plot( x_locs_lines[i], np.array(powers_lp_lines[i])/np.array(powers_lines[i]), 
             ls=lstyle[i], linewidth=4, color=colors[i] )
    if addmarker:
        add_sep_pfmax(i, x_locs_lines[i], np.array(powers_lp_lines[i])/np.array(powers_lines[i]))
    
plot_style()
ax.set_ylabel(r"$\int \mathcal{P}(f) \mathrm{d} f_{|St<"+f"{cutoff_st}"+r"}/\int \mathcal{P}(f) \mathrm{d} f$" )
ax.set_ylim([0,1])
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
figname = "pf_power_lp_ratio"
plt.savefig( figname + '.png' )
plt.close()

# =============================================================================
# plot the ratio of low-frequency-pass power spectra to the high-frequency-pass power spectra

fig, ax = plt.subplots( figsize=(15, 8), constrained_layout=True )
for i in range( len(cases) ):
    ax.plot( x_locs_lines[i], np.array(powers_lp_lines[i])/np.array(powers_hp_lines[i]), 
             ls=lstyle[i], linewidth=4, color=colors[i] )
    if addmarker:
        add_sep_pfmax(i, x_locs_lines[i], np.array(powers_lp_lines[i])/np.array(powers_hp_lines[i]))
    
plot_style()
ax.set_ylabel(r"$\int \mathcal{P}(f) \mathrm{d} f_{|St<"+f"{cutoff_st}"+r"}/\int \mathcal{P}(f) \mathrm{d} f_{|St>"+f"{cutoff_st}"+r"}$" )
ax.set_ylim([0,3.1])
ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
figname = "pf_power_l2h_ratio"
plt.savefig( figname + '.png' )
plt.close()