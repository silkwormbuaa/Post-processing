# -*- coding: utf-8 -*-
'''
@File    :   plt_style.py
@Time    :   2023/05/09 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import pickle

import numpy             as     np

import pandas            as     pd

import matplotlib.pyplot as     plt

import matplotlib.ticker as     ticker

from   matplotlib.patches import Circle

from   mpl_toolkits.axes_grid1 import make_axes_locatable


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
# >>> Plot 2d dmd mode                                            (Nr.)
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
                   St=None,
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
    
    # print St one figures
    
    if St is not None:
        plt.text(2,3, f'St={St:.3f}', fontsize=18)
    
    # show colorbar?
    
    if colorbar is True:
        cbar = plt.colorbar(cs, cax=ax.inset_axes([1.1,0.3,0.03,0.5]))
        cbar.ax.set_ylabel(r'$\Re({\phi}_p)$',fontsize=15)
    
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
# >>> plot combined 2d dmd mode                                   (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/26  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_combined_dmd_mode( grids1, v1, dir1,
                            grids2, v2, dir2,
                            figsize=None,
                            filename=None,
                            pure=None,
                            colorbar=False,
                            clevel=None,
                            title=None):
    
    if figsize is None:
        
        # set default figure size
        figsize = (30,20)
  
    fig = plt.figure(figsize=figsize)
    ax  = fig.add_subplot( 111, projection='3d' )
    
    # level of color bar
    
    if clevel is None:
        
        clevel = np.linspace( min( np.min(v1), np.min(v2) ),
                              max( np.max(v1), np.max(v2) ),
                              51 )
           
    # plane direction to determine passing value
    # since ax.view_init() only support two axes to adjust view,
    # in order to show a better view angle, switch y,z values.    
      
    if dir1 == 'x':
        X = v1;        Y = grids1[0]; Z = grids1[1]
    elif dir1 == 'y':
        X = grids1[0]; Y = grids1[1]; Z = v1 ; dir1 = 'z'
    elif dir1 == 'z':
        X = grids1[0]; Y = v1; Z = grids1[1] ; dir1 = 'y'
    
    # intersection line between Y,Z planes
        
    line = np.array([[np.min(X),0,0],[np.max(X),0,0]])
    
    ax.plot(line[:,0],line[:,1],line[:,2],
            color='black',linewidth=1.0,zorder=100)
        
    contourf1 = ax.contourf( X,Y,Z, 
                             zdir=dir1, 
                             alpha=1.0,
                             levels=clevel, 
                             offset=0,
                             cmap='bwr' )

    if dir2 == 'x':
        X = v2;        Y = grids2[0]; Z = grids2[1]
    elif dir2 == 'y':
        X = grids2[0]; Y = grids2[1]; Z = v2; dir2='z'
    elif dir2 == 'z':
        X = grids2[0]; Y = v2; Z = grids2[1]; dir2='y'
        
    ax.contourf( X,Y,Z, 
                 zdir = dir2, 
                 alpha = 0.9,
                 offset = 0,
                 levels=clevel, 
                 cmap='bwr' )
    
    # view angle    
    ax.view_init( elev=30.0, azim=-120.0 )
    
    # parameters for axes
    
    ax.set_xlabel(r'$(x-x_{imp})/{\delta}_0$',fontsize=40)
    ax.set_zlabel(r'$y_s/{\delta}_0$',fontsize=40)
    ax.set_ylabel(r'$z/{\delta}_0$',fontsize=40)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    ax.zaxis.set_major_locator(ticker.MultipleLocator(2.0))
    
    ax.tick_params(axis='x',labelsize=30,pad=20)
    ax.tick_params(axis='y',labelsize=30,pad=10)
    ax.tick_params(axis='z',labelsize=30,pad=10)

    ax.xaxis.labelpad=70
    ax.yaxis.labelpad=20
    ax.zaxis.labelpad=20

    ax.set_zlim(0,7)
    ax.set_ylim(-2,2)
    
    # colorbar
    
    if colorbar is True: 
        cbar = plt.colorbar(contourf1, 
                            cax=ax.inset_axes([1.05,0.3,0.03,0.5]))
        cbar.ax.set_ylabel(r'$\Re({\phi}_p)$',fontsize=40)
        cbar.ax.tick_params(labelsize=30)

    # title        
    if title is not None: 
        ax.text(0.0,8.0,5.0, title, fontsize=40, ha='center')
        
    # set aspect of bounding box
    ax.set_box_aspect([150,40,50])

    plt.tight_layout()  # tried to remove rounding blank

    if filename:
        
        plt.savefig( filename, bbox_inches='tight' ) # tight bouding box
        print(f"{filename} is output.\n")
    
    plt.close()
    


# ----------------------------------------------------------------------
# >>> Plot Z slice                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/12  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_slicez_stat( xx, yy, v, 
                      sonic=True, 
                      separation=True,
                      figsize=None, 
                      filename=None,
                      col_map=None,
                      cbar_label=None,
                      cbar_levels=None):
    
    if figsize is None:
        figsize = (15,8)
    
    fig, ax = plt.subplots( figsize=figsize )
    
    if col_map is None: col_map='viridis'
    
    if cbar_levels is None: cbar_levels=51
    
    cs = ax.contourf( xx, yy, v, cmap=col_map, levels=cbar_levels )
    
    if sonic: 
        with open('soniclines.pkl','rb') as f:
            lines = pickle.load( f )
        
        for line in lines:
            x_sonic = line[:,0]
            y_sonic = line[:,1]
            ax.plot(x_sonic,y_sonic,'white',linewidth=0.8)


    if separation:
        with open('separationlines.pkl','rb') as f:
            lines = pickle.load( f )
            
        for line in lines:
            x_sep = line[:,0]
            y_sep = line[:,1]
            ax.plot(x_sep,y_sep,'red',linewidth=0.8)

    cbar = plt.colorbar(cs)
    cbar.ax.set_ylabel(cbar_label,fontsize=15)
    cbar.ax.tick_params(labelsize=15)

    ax.tick_params(axis='x',labelsize=15,pad=10)
    ax.tick_params(axis='y',labelsize=15,pad=10)


    # set axises stay with contour and x,y unit length equal
    
    plt.axis('tight')
    plt.gca().set_aspect('equal', adjustable='box')

    if filename:
        plt.savefig( filename )
        print(f"{filename} is output.\n")

    plt.close()        


# ----------------------------------------------------------------------
# >>> Plot X slice                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/12  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_slicex_stat( zz, yy, v,
                      wall=True, 
                      sonic=True, 
                      separation=False,
                      figsize=None, 
                      filename=None,
                      col_map=None,
                      cbar_label=None):
    
    if figsize is None:
        figsize = (15,15)
    
    fig, ax = plt.subplots( figsize=figsize )
    
    if col_map is None: col_map='viridis'
    
    cs = ax.contourf( zz, yy, v, cmap=col_map, levels=51 )
    
    if wall:
        df_wall = pd.read_csv( 'wall_X.dat', delimiter=r'\s+' )
        z_wall = np.array(df_wall['x'])
        y_wall = np.array(df_wall['y'])
#        ax.plot(z_wall,y_wall,'black')
        ax.fill_between(z_wall,-0.1,y_wall,color='grey')
    
    
    if sonic: 
        with open('soniclines.pkl','rb') as f:
            lines = pickle.load( f )
        
        for line in lines:
            x_sonic = line[:,0]
            y_sonic = line[:,1]
            ax.plot(x_sonic,y_sonic,'white',linewidth=0.8)
    
    
    if separation:
        with open('separationlines.pkl','rb') as f:
            lines = pickle.load( f )
            
        for line in lines:
            x_sep = line[:,0]
            y_sep = line[:,1]
            ax.plot(x_sep,y_sep,'red',linewidth=0.8)

    cbar = plt.colorbar(cs)
    cbar.ax.set_ylabel(cbar_label,fontsize=15)
    cbar.ax.tick_params(labelsize=15)

    ax.tick_params(axis='x',labelsize=15,pad=10)
    ax.tick_params(axis='y',labelsize=15,pad=10)


    # set axises stay with contour and x,y unit length equal
    
    plt.axis('tight')
    plt.gca().set_aspect('equal', adjustable='box')

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