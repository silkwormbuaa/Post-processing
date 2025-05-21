#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   mimic_shock.py
@Time    :   2025/04/28 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   mimic the shock motion and check the influence on statistics
'''

import os
import sys
import numpy                as     np
import matplotlib.pyplot    as     plt
from   matplotlib.animation import FuncAnimation

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.directories    import create_folder


def main():

    outfolder = '/home/wencan/temp/DataPost/midRe/mimic_shock/delta_compare'
    os.chdir( create_folder( outfolder ) )
    save = True
    
    t = 0.0
    amplitude = 1.0
    delta     = 0.3
    delta2    = 0.2
    x_c       = 0.0
    x_amp     = 0.3
    
    x, p = pressure_shape( t, amplitude, delta, x_c, x_amp )
    plt.plot( x, p )
    plt.xlabel( r'$x$' )
    plt.ylabel( r'$p$' )
    plt.title( r'$t=0.0$' )
    plt.grid()
    if save:
        plt.savefig( 'pressure_t0.png' )
    else: plt.show()
    plt.close()
    
    
    animate_pressure(t, amplitude, delta, x_c, x_amp, save)
    
    
    pressure_array  = list()
    pressure_array2 = list()
    
    for t in np.linspace(0.0, 1.0, 1000, endpoint=False):
        x, p = pressure_shape( t, amplitude, delta, x_c, x_amp )
        pressure_array.append( p )

    for t in np.linspace(0.0, 1.0, 1000, endpoint=False):
        x, p = pressure_shape( t, amplitude, delta2, x_c, x_amp )
        pressure_array2.append( p )
    
    pressure_array = np.array(pressure_array)
    pressure_array = np.transpose(pressure_array)
    mean_pressure  = np.mean(pressure_array, axis=1)
    pressure_fluc  = pressure_array - mean_pressure[:, np.newaxis]
    rms_pressure   = np.std(pressure_fluc, axis=1)
    pressure_grad  = np.gradient(mean_pressure, x)

    pressure_array2 = np.array(pressure_array2)
    pressure_array2 = np.transpose(pressure_array2)
    mean_pressure2  = np.mean(pressure_array2, axis=1)
    pressure_fluc2  = pressure_array2 - mean_pressure2[:, np.newaxis]
    rms_pressure2   = np.std(pressure_fluc2, axis=1)
    pressure_grad2  = np.gradient(mean_pressure2, x)    

    # plt.plot(x, mean_pressure)
    # plt.plot(x, mean_pressure2)
    # plt.xlabel( r'$x$' )
    # plt.ylabel( r'$\langle p \rangle$' )
    # plt.grid()
    # if save: plt.savefig( 'mean_pressure.png' )
    # else:    plt.show()
    # plt.close()

    # plt.plot(x, rms_pressure)
    # plt.plot(x, rms_pressure2)
    # plt.xlabel( r'$x$' )
    # plt.ylabel( r'$p^{\prime}_{r.m.s.}$' )
    # plt.grid()
    # if save: plt.savefig( 'rms_pressure.png' )
    # else:    plt.show()
    # plt.close()

    # plt.plot(x, pressure_grad)
    # plt.plot(x, pressure_grad2)
    # plt.xlabel( r'$x$' )
    # plt.ylabel( r'$\frac{\partial \langle p \rangle}{\partial x}$' )
    # plt.grid()
    # if save: plt.savefig( 'pressure_grad.png' )
    # else:    plt.show()
    # plt.close()
    
    # extract the pressure at a location before the shock mean location
    
    pres_up = pressure_array2[700,:]
    
    # add some noise to the pressure
    pres_up += np.random.normal(0.0, 0.02, size=pres_up.shape)
    
    pres_up = (pres_up - np.mean(pres_up)) / np.std(pres_up)
    
    plt.plot(np.linspace(0.0, 1.0, 1000), pres_up)
    plt.hlines( np.mean(pres_up), 0.0, 1.0, color='black', linestyle='--', label='mean')
    plt.xlabel( r'$t$' )
    plt.ylabel( r'$p$' )
    plt.grid()
    plt.show()
    plt.close()
    
    # plot the pres_up's p.d.f 
    # normalize the pres_up first
    print(f"mean: {np.mean(pres_up)}, std: {np.std(pres_up)}")
    hist, bin_edges = np.histogram( pres_up, bins=50, range=(-3,3 ), density=True )
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:]) 
    
    # - normal distribution
    
    x = np.linspace(-3,3,100)
    y = 1/np.sqrt(2*np.pi) * np.exp(-0.5*x**2)
    
    
    fig, ax = plt.subplots( figsize=(8,6) )
    
    ax.plot( bin_centers, hist, 'o', markersize=10 )
    ax.plot( x, y, 'gray', ls=':', lw=1.5 )
    
    plt.xlabel( r'$\sigma_{pw}$' )
    plt.ylabel( r'$p.d.f.$' )
    plt.show()
    plt.close()
    
    
def animate_pressure(t, amplitude, delta, x_c, x_amp, save=False):
    
    fig, ax = plt.subplots( figsize=(12, 6) )
    
    x, p = pressure_shape( t, amplitude, delta, x_c, x_amp )
    line, = ax.plot( x, p )
    ax.grid()
    
    def update(frame):
        nonlocal t
        t += 0.01
        x, p = pressure_shape( t, amplitude, delta, x_c, x_amp )
        line.set_ydata(p)
        return line,
    
    ani = FuncAnimation(fig, update, frames=100, interval=50)
    if save:
        ani.save('shock_wave.mp4', fps=25, extra_args=['-vcodec', 'libx264'])
    else:
        plt.show()
    plt.close()

def pressure_shape( t, amplitude, delta, x_c, x_amp ):
    
    """
    t:         the time
    amplitude: the pressure jump amplitude
    delta:     the width of the pressure jump
    x_c:       the center of the pressure jump
    x_amp:     the streamwise meandering amplitude
    """
    
    p0 = 1.0
    
    x = np.linspace( -1.0, 1.0, 1000 )
    
    x0 = x_c + x_amp * np.sin( 2.0*np.pi*t )   # freq = 1 Hz
    
    p = p0 + 0.5 * amplitude * (1 + np.tanh( (x-x0)/delta ) )

    return x, p
    

def plot_pressure( t, x, p ):
    
    fig, ax = plt.subplots( figsize=(12, 6) )
    
    ax.plot( x, p, label=f't={t:.2f}s' )
    
    ax.set_xlabel( r'$x$' )
    ax.set_ylabel( r'$p$' )
    ax.legend()
    ax.grid()
    
    plt.show()
    plt.close()




# =============================================================================
if __name__ == "__main__":

    main()

