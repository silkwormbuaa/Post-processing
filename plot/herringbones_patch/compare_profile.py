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

import matplotlib.pyplot  as     plt
import matplotlib.ticker  as     ticker
import matplotlib.markers as     markers
import numpy              as     np
import pandas             as     pd

source_dir = os.path.realpath(__file__).split('plot')[0]
sys.path.append( source_dir )

from   vista.directories  import Directories
from   vista.line         import ProfileData

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 20

def main():
    loc    = -12  # -12
    outfolder = '/home/wencan/temp/DataPost/herringbones_patch/profile_-12'
    cases  = ['smooth_adiabatic','250821','250710']
    colors = ['black'           ,'red'   ,'blue'  ]
    
    vars   = ['u','v','T','rho']
    labels = [r"$\langle u \rangle /u_{\infty}$",
              r"$\langle v \rangle /u_{\infty}$",
              r"$\langle T \rangle /T_{0}$",
              r"$\langle \rho \rangle /\rho_{\infty}$"]
    
    ranges = [[0.0,1.0],[-0.05,0.05],[0.5,1.0],[0.4,1.0]]
    norms  = [507,507,288.2,0.9886]
    
    casepaths = [ '/home/wencan/temp/'+case for case in cases] 
    
    linems, lineds, linecs = read_case_profiles( casepaths[0], loc )
    linem1, lined1, linec1 = read_case_profiles( casepaths[1], loc )
    linem2, lined2, linec2 = read_case_profiles( casepaths[2], loc )
    
    os.chdir( outfolder )
    for i, var in enumerate( vars ): 
    
        fig, ax = plt.subplots( figsize=[6,6] )
        fig.subplots_adjust( left=0.18, right=0.95, top=0.95, bottom=0.20 )
        ax.plot(linems.df[var]/norms[i],linems.df['y']/5.2,colors[0],'-')
        ax.plot(linem1.df[var]/norms[i],linem1.df['y']/5.2,colors[1],'-')
        ax.plot(linem2.df[var]/norms[i],linem2.df['y']/5.2,colors[2],'-')

        ax.plot(lined1.df[var]/norms[i],lined1.df['y']/5.2,colors[1],ls='--')
        ax.plot(lined2.df[var]/norms[i],lined2.df['y']/5.2,colors[2],ls='--')
        
        ax.plot(linec1.df[var]/norms[i],linec1.df['y']/5.2,colors[1],ls=':')
        ax.plot(linec2.df[var]/norms[i],linec2.df['y']/5.2,colors[2],ls=':')
        
        ax.minorticks_on()
        ax.set_xlabel( labels[i] )
        ax.set_xlim( ranges[i] )
        
        ax.set_ylim([0.0,1.2])
        yticks = np.linspace(0, 1.0, 6)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f'{y:.1f}' for y in yticks])
        ax.set_ylabel( r'$y/\delta_0$' )


        ax.xaxis.set_major_locator( ticker.MaxNLocator(6) )
        ax.xaxis.set_minor_locator( ticker.AutoMinorLocator(5) )
        ax.yaxis.set_major_locator( ticker.MaxNLocator(6) )
        ax.yaxis.set_minor_locator( ticker.AutoMinorLocator(5) )

        ax.tick_params(which='major',
                    axis='both',
                    direction='out',
                    length=15,
                    width=1.5)
        ax.tick_params(which='minor',
                    axis='both', 
                    direction='out',
                    length=10,
                    width=1.0)
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(2)
        
        plt.savefig( f'profile_{var}.png', dpi=300 )
        print(f"Figure saved to {os.getcwd()}/profile_{var}.png.")
#        plt.show()
        plt.close()

def read_case_profiles( casepath, loc ):
    
    def read_line( file ):
        line = ProfileData( file )
        line.shift_y( 0.0 )
        line.inner_normalize( 'wall_statistics.dat' )
        line.vd_transform()

        return line
    
    dirs = Directories( casepath )
    datafolder = dirs.pp_statistics + f'/profile_{str(int(loc))}'
    os.chdir( datafolder )

    linem = read_line(  'profile_mean.dat'   )
    lined = read_line( f'profile_{str(int(loc))}_DL.dat' )
    linec = read_line( f'profile_{str(int(loc))}_CL.dat' )
    
    return linem, lined, linec

# =============================================================================
if __name__ == "__main__":

    main()

