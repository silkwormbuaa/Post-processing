# -*- coding: utf-8 -*-
'''
@File    :   plt_style.py
@Time    :   2023/05/09 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import numpy             as     np

import pandas            as     pd

import matplotlib.pyplot as     plt

import matplotlib.ticker as     ticker

from   matplotlib.patches import Circle


# ----------------------------------------------------------------------
# >>> Plot eigen values around a unit circle                     (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_eigens( eigens,
                 eigens2=None,
                 figsize=None, 
                 filename=None, 
                 show_circle=True,
                 set_view=False,
                 pure=False):
    
    if eigens is None:
        
        raise ValueError('Please compute eigen values first.')
    
    
    if figsize is None:
        
        # Set default figure size
        
        figsize = (10,10) 
        
        
    
    fig, ax = plt.subplots( figsize=figsize )
    
    # Show unit circle
    
    if show_circle:
        
        circle = Circle( (0,0), 
                         radius=1, 
                         color='green', 
                         fill=False, 
                         ls='--',
                         label='unit circle')

        ax.add_artist( circle )
        
        
    # Plot the first eigen values
    
    ax.scatter( eigens.real, 
                eigens.imag, 
                marker='o',
                edgecolors='gray',
                facecolors='none',
                label='Eigenvalues')
    
    
    # Plot the compared eigen values set if they are given
    
    if eigens2 is not None:
        
        ax.scatter( eigens2.real, 
                    eigens2.imag, 
                    marker='x',
                    c='blue',
                    label='Eigenvalues2')
    
    # Set_view
    
    if set_view:
        
        ax.set_xlim([ -1.2, 1.2 ])
        ax.set_ylim([ -1.2, 1.2 ])
     
    # Ticks

    ax.minorticks_on()
    
    ax.tick_params( which='major',
                    axis='both',
                    direction='in',
                    length=10,
                    width=2)
    
    ax.tick_params( which='minor',
                    axis='both', 
                    direction='in',
                    length=5,
                    width=1)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    
    ax.tick_params(axis='x',labelsize=15)
    ax.tick_params(axis='y',labelsize=15)
    
    # Grids
    
    ax.grid(which='major', ls=':', linewidth=1)
    ax.grid(which='minor', ls=':', linewidth=0.5)
    
    # Labels
    
    ax.set_xlabel( "Real part",fontdict={'size':20})
    ax.set_ylabel( "Imaginary part",fontdict={'size':20})
    
    # With default labels or get pure figure
    
    if pure: fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
    # Save figure
    
    if filename: 
        
        plt.savefig( filename )
        print(f"{filename} is saved.\n")
    


# ----------------------------------------------------------------------
# >>> Plot amplitude_St                                          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_amp_st( st, amp1, amp2=None, 
                 figsize=None, 
                 filename=None, 
                 pure=False, 
                 hidesmall=True ):
    
    if amp1 is None:
        
        raise ValueError('Please compute amplitude first.')
    
    
    if figsize is None:
        
        # Set default figure size
        
        figsize = (10,10) 
        
    
    fig, ax = plt.subplots( figsize=figsize )      
        
    # Plot the first amplitudes
    
    ax.scatter( st, 
                amp1, 
                marker='o',
                edgecolors='gray',
                facecolors='none',
                label='Standard' )
    
    if amp2 is not None:
        
        # if hide small, drop small amplitude less than 1
        
        df = pd.DataFrame( st, columns=['st'])
        
        df['amp2'] = amp2
            
        if hidesmall:
            
            df = df.drop( df[ abs(df['amp2'])<0.001 ].index )
        
        
        ax.scatter( np.array( df['st'] ), 
                    np.array( df['amp2'] ), 
                    marker='x',
                    c='red',
                    label='Sparsity-promoting' )
    
     
    # Ticks

    ax.minorticks_on()
    
    ax.tick_params( which='major',
                    axis='both',
                    direction='in',
                    length=10,
                    width=2)
    
    ax.tick_params( which='minor',
                    axis='both', 
                    direction='in',
                    length=5,
                    width=1)
    
    ax.set_yscale( "symlog", linthresh = 1 )
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5.0))
#    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    
    ax.tick_params(axis='x',labelsize=15)
    ax.tick_params(axis='y',labelsize=15)
    
    # Grids
    
    ax.grid(which='major', ls=':', linewidth=1)
    ax.grid(which='minor', ls=':', linewidth=0.5)
    
    # Labels
    
    ax.set_xlabel( r"$St_{L_{sep}}$",fontdict={'size':20})
    ax.set_ylabel( "amplitude",fontdict={'size':20})
    
    # With default labels or get pure figure
    
    if pure: fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
    # Save figure
    
    if filename: 
        
        plt.savefig( filename )
        print(f"{filename} is saved.\n")

    plt.close() # must be added, otherwise will give too many figures warning.


# ----------------------------------------------------------------------
# >>> Plot psi(relative amp)_St                                       (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_psi_st( st, psi1, psi2=None, 
                 figsize=None, 
                 filename=None, 
                 pure=False, 
                 hidesmall=True,
                 xlim=None,
                 gray=None):
    
    if psi1 is None:
        
        raise ValueError('Please compute amplitude first.')
    
    
    if figsize is None:
        
        # Set default figure size
        
        figsize = (10,10) 
        
    
    fig, ax = plt.subplots( figsize=figsize )      
        
    # Plot the first amplitudes
    df = pd.DataFrame( st, columns=['st'])
    
    df['psi1'] = psi1
    
    if gray is not None: df['gray'] = gray
    
    if psi2 is not None: df['psi2'] = psi2
    
    df = df.drop( df[ df['st']<0.0 ].index )
    

    # Ticks

#    ax.minorticks_on()
    
    ax.tick_params( which='major',
                    axis='both',
                    direction='out',
                    length=10,
                    width=2)
    
    ax.tick_params( which='minor',
                    axis='both', 
                    direction='out',
                    length=5,
                    width=1)
    
    
#    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_tick_params( which='both', zorder=-1 )
    ax.tick_params( axis='x', labelsize=15 )
    ax.tick_params( axis='y', labelsize=15 )
    
    # axis limit
    
    if xlim:
        ax.set_xscale( "log" )
        ax.set_xlim( xlim )
    
    else: ax.set_xlim( [0.0,3.0] )
    
    ax.set_yscale( "log" )
    ax.set_ylim( [0.0001,1.0] )
    
    # Grids
    
    ax.grid(which='major', ls=':', linewidth=1)
    ax.grid(which='minor', ls=':', linewidth=0.5)
    
    # Labels
    
    ax.set_xlabel( r"$St_{L_{sep}}$",fontdict={'size':20})
    ax.set_ylabel( r"$|\psi_{k}|$",fontdict={'size':20})
    
    # Plot line
    
    if gray is not None:
        
        st1 = np.array(df['st'])
        v_psi1 = np.array(df['psi1'])
        
        for i, value in enumerate( df['gray'] ):
            plt.vlines( st1[i],
                        0.0,
                        v_psi1[i],
                        linewidth=1,
                        alpha= value,
                        colors='black')
    
    else:
    
        plt.vlines( df['st'],
                    0.0,
                    df['psi1'], 
                    linewidth=1,
                    alpha=1.0,
                    colors='black')
    
    if psi2 is not None:
        
        # if hide small, drop small amplitude less than 1
            
        if hidesmall:
            
            df2 = df.drop( df[ abs(df['psi2'])<1.0e-8 ].index )
        
        
        plt.vlines( np.array( df2['st'] ),
                    0.0, 
                    np.array( df2['psi2'] ), 
                    linewidth=1.5,
                    colors='red',
                    alpha=1.0,
                    label='Sparsity-promoting' )
        
        plt.plot( np.array( df2['st']),
                  np.array( df2['psi2']),
                  'o',
                  color='red',
                  markersize=4)
    # Save figure
    
    if filename: 
        
        plt.savefig( filename )
        print(f"{filename} is saved.\n")
        
    # With default labels or get pure figure
    
    if pure: 
        
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        plt.savefig( 'apure'+filename )
        print(f"{'apure'+filename} is saved.\n")
    
    plt.close()
        


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/06/09  - created
#
# Desc
#   - 2D contours of dmd modes
# ----------------------------------------------------------------------

def plot_dmd_mode( grids, 
                   v,
                   figsize=None, 
                   filename=None,
                   pure=None,
                   colorbar=False,
                   clevel=None,
                   title=None):
    
    if (grids is None) or (v is None):
        
        raise ValueError('grids or mode value is None!')
    
    if figsize is None:
        
        # Set default figure size
        
        figsize = (15,4)
        
    
    fig, ax = plt.subplots( figsize = figsize )
    
    if clevel is None:
        
        clevel = np.linspace( np.min(v), np.max(v), 51 )
    
    
    cs = ax.contourf(grids[0],
                     grids[1],
                     v,
                     levels=clevel,
                     cmap='bwr')
    
    # show colorbar?
    
    if colorbar is True: plt.colorbar( cs )
    
    # set title
    
    if title: ax.set_title(f"{title}",size=20)
    
    # Set view
    
    ax.set_xlim([np.min(grids[0]),np.max(grids[0])])
    ax.set_ylim([np.min(grids[1]),np.max(grids[1])])
    
    # set axises stay with contour and x,y unit length equal
    
    plt.axis('tight')
    plt.gca().set_aspect('equal', adjustable='box')
    
    # save figures
    
    if filename:
        
        plt.savefig( filename )
        
        print(f"{filename} is output.\n")
    
    plt.close()



# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

 
    np.random.seed(123)
    real_part = np.random.random(10)*2-1
    imag_part = np.random.random(10)*2-1
    eigenvalues = real_part + 1j * imag_part
    
    np.random.seed(111)
    real_part = np.random.random(10)*2-1
    imag_part = np.random.random(10)*2-1
    eigenvalues2 = real_part + 1j * imag_part

    
    plot_eigens(eigenvalues,eigenvalues2,set_view=True)



# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()