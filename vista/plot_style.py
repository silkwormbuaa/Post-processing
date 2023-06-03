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
    
    if filename: plt.savefig( filename )
    
    
    plt.show()



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
            
            df = df.drop( df[ df['amp2']<0.01 ].index )
        
        
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
    
    if filename: plt.savefig( filename )
    
    
    plt.show()


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