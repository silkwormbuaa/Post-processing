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

import matplotlib
import matplotlib.pyplot as     plt
import matplotlib.colors as     colors
import matplotlib.ticker as     ticker

from   matplotlib.patches import Circle

from   mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size']   = 40

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
                 xlim=None,
                 ylim=None):
    
    if eigens is None:
        
        raise ValueError('Please compute eigen values first.')
    
    
    if figsize is None:
        
        # Set default figure size
        
        figsize = (7,5.5) 
        
    
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
                label='Eigenvalues',
                s=40)
    
    
    # Plot the compared eigen values set if they are given
    
    if eigens2 is not None:
        
        ax.scatter( eigens2.real, 
                    eigens2.imag, 
                    marker='x',
                    c='red',
                    label='Eigenvalues2',
                    s=150,
                    linewidth=2)
    
    # Set_view
    
    if set_view:
        
        ax.set_xlim([ -1.2, 1.2 ])
        ax.set_ylim([ -1.2, 1.2 ])
        
    if xlim is not None: ax.set_xlim( xlim )
    if ylim is not None: ax.set_ylim( ylim )
     
    # Ticks

    ax.minorticks_on()
    
    ax.tick_params( which='major',
                    axis='both',
                    direction='in',
                    length=20,
                    width=2)
    
    ax.tick_params( which='minor',
                    axis='both', 
                    direction='in',
                    length=10,
                    width=1.5)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.04))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')
    
    # Grids
    
#    ax.grid(which='major', ls=':', linewidth=1)
#    ax.grid(which='minor', ls=':', linewidth=0.5)
    
    # Labels
    
    ax.set_xlabel( r"$\mathfrak{R}(\mu_i)$")
    ax.set_ylabel( r"$\mathfrak{I}(\mu_i)$")
    
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
    # adjust the content region
    
    fig.subplots_adjust(left=0.32, right=0.95, bottom=0.25, top=0.95)
    
    # set background transparent
    
    fig.patch.set_alpha(0.0)
    
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
        
        figsize = (10,8) 
        
    
    fig, ax = plt.subplots( figsize=figsize,constrained_layout=True )      
        
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
    
    # Grids
    
#    ax.grid(which='major', ls=':', linewidth=1)
#    ax.grid(which='minor', ls=':', linewidth=0.5)
    
    # Labels
    
    ax.set_xlabel( r"$St_{L_{sep}}$")
    ax.set_ylabel( "amplitude")
    
    # With default labels or get pure figure
    
    if pure: fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)
    
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
                 hidesmall=True,
                 xlim=None,
                 gray=None):
    
    if psi1 is None:
        
        raise ValueError('Please compute amplitude first.')
    
    
    if figsize is None:
        
        # Set default figure size
        
        figsize = (15,9) 
        
    fig, ax = plt.subplots( figsize=figsize )      
        
    # Plot the first amplitudes
    df = pd.DataFrame( st, columns=['st'])
    
    df['psi1'] = psi1
    
    if gray is not None: df['gray'] = gray
    
    if psi2 is not None: df['psi2'] = psi2
    
    df = df.drop( df[ df['st']<0.0 ].index )
    

    # Ticks

    ax.minorticks_on()
    
    ax.tick_params( which='major',
                    axis='both',
                    direction='out',
                    length=20,
                    width=2)
    
    ax.tick_params( which='minor',
                    axis='both', 
                    direction='out',
                    length=10,
                    width=1.5)

    # set major tick on y axis is 0.0001,0.001,0.01,0.1,1.0
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=5))

#    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_tick_params( which='both', zorder=-1 )
    
    # axis limit
    
    if xlim:
        ax.set_xscale( "log" )
        ax.set_xlim( xlim )
    
    else: ax.set_xlim( [0.0,3.0] )
    
    ax.set_yscale( "log" )
    ax.set_ylim( [0.0001,1.0] )
    
    # Grids
    
#    ax.grid(which='major', ls=':', linewidth=1)
#    ax.grid(which='minor', ls=':', linewidth=0.5)
    
    # Labels
    
    ax.set_xlabel( r"$St_{L_{sep}}$")
    ax.set_ylabel( r"$|\psi_{i}|$")
    
    # hide x axis label and tick number
    ax.xaxis.set_tick_params( which='both', labelbottom=True )
    
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
                    linewidth=3,
                    colors='red',
                    alpha=1.0,
                    label='Sparsity-promoting' )
        
        plt.plot( np.array( df2['st']),
                  np.array( df2['psi2']),
                  'o',
                  color='red',
                  markersize=8)
        
    # set the bounding box of axes
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(3)

    # Save figure
    
    if filename: 
        
        fig.subplots_adjust(left=0.125, right=0.9375, bottom=0.1667, top=0.9167)
        plt.savefig( filename )
        print(f"{filename} is saved.\n")
        

    
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
                            shock_shape=None,
                            colorbar=False,
                            clevel=None,
                            title=None):
    
    if figsize is None:
        
        # set default figure size
        figsize = (20,15)
  
    fig = plt.figure(figsize=figsize)
    ax  = fig.add_subplot( 111, projection='3d' )
    
    # level of color bar
    
    if clevel is None:
        
        clevel = np.linspace( min( np.min(v1), np.min(v2) ),
                              max( np.max(v1), np.max(v2) ),
                              36 )
           
    # plane direction to determine passing value
    # since ax.view_init() only support two axes to adjust view,
    # in order to show a better view angle, switch y,z values.    
      
    if dir1 == 'x':
        X = v1;        Y = grids1[0]; Z = grids1[1]
    elif dir1 == 'y':
        X = grids1[0]; Y = v1;        Z = grids1[1]
    elif dir1 == 'z':
        X = grids1[0]; Y = grids1[1]; Z = v1
    
    # intersection line between Y,Z planes
        
    line = np.array([[np.min(X),0,0],[np.max(X),0,0]])
    
    ax.plot(line[:,0],line[:,1],line[:,2],
            color='black',linewidth=1.0,zorder=100)
        
    contourf1 = ax.contourf( X,Y,Z, 
                             zdir=dir1, 
                             alpha=1.0,
                             levels=clevel, 
                             offset=0,
                             cmap='RdBu_r',
                             extend='both')
    
    if shock_shape is not None:
        lines = pickle.load( open(shock_shape,'rb') )
        for line in lines:
            x_shock = line[:,0]
            y_shock = line[:,1]
            ax.plot(x_shock,np.zeros_like(x_shock),y_shock,zdir=dir1,
                    color='black',linewidth=1.0,zorder=101)

    if dir2 == 'x':
        X = v2;        Y = grids2[0]; Z = grids2[1]
    elif dir2 == 'y':
        X = grids2[0]; Y = v2;        Z = grids2[1]
    elif dir2 == 'z':
        X = grids2[0]; Y = grids2[1]; Z = v2
        
    ax.contourf( X,Y,Z, 
                 zdir = dir2, 
                 alpha = 0.9,
                 offset = 0,
                 levels=clevel, 
                 cmap='RdBu_r',
                 extend='both')
    
    # view angle    
    ax.view_init( elev=25, azim=-32, vertical_axis='y' )
    ax.set_proj_type('ortho')
    
    # parameters for axes
    
    ax.set_xlabel(r'$(x-x_{imp})/{\delta}_0$', labelpad=100)
    ax.set_ylabel(r'$y_s/{\delta}_0$', labelpad=20 )
    ax.set_zlabel(r'$z/{\delta}_0$', labelpad=30 )
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    ax.zaxis.set_major_locator(ticker.MultipleLocator(2.0))
    
    ax.tick_params( which='major',
                    axis='x',
                    direction='out',
                    length=15,
                    width=2)
    
    ax.tick_params(axis='x',pad=10)
    ax.tick_params(axis='y',pad=5)
    ax.tick_params(axis='z',pad=8)

    ax.set_zlim(-2,2)
    ax.set_ylim(0,8)
    
    # add 3D axes to show the view angle
    
    ax.quiver(-20, -5, 0, 1, 0, 0, color='r', label='X-axis')
    ax.quiver(-20, -5, 0, 0, 1, 0, color='g', label='Y-axis')
    ax.quiver(-20, -5, 0, 0, 0, 1, color='b', label='Z-axis')
    
    # colorbar
    
    if colorbar is True: 
        cbar = plt.colorbar(contourf1, 
 #                           cax=ax.inset_axes([0.0,0.3,0.03,0.3]),
                            orientation='vertical', 
                            shrink=0.3,
                            location='left',
                            aspect=10)
        
        cbar.outline.set_linewidth(2)
        
        cbar.ax.set_xlabel( r"$\mathfrak{R}(\phi_p)$",
                            fontsize=40,
                            labelpad=25,
                            loc='center')
        
        cbar.ax.tick_params( direction='in',
                             left=True,right=False,
                             labelleft=True,labelright=False,
                             length=20.0,
                             width=2.0)

    # title        
    if title is not None: 
        ax.text(0.0,15.0,0.0, title, ha='center')
        
    # set aspect of bounding box
    ax.set_box_aspect([6,50,12])
    fig.patch.set_facecolor('white')

#    plt.tight_layout()  # tried to remove rounding blank
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.1, top=0.9)
    fig.tight_layout()

    if filename:
        
        plt.savefig( filename ) # tight bouding box
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
                      sonic=None, 
                      separation=None,
                      boundary=None,
                      shockshape=None,
                      DS=None,
                      figsize=None, 
                      filename=None,
                      col_map=None,
                      cbar_label=None,
                      cbar_levels=None,
                      cbar_ticks=None,
                      tag=None,
                      tag_loc=None,
                      x_lim=None,
                      y_lim=None,
                      wall=None,
                      pure=False):
    
    if figsize is None:
        figsize = (20,10)
    
    fig, ax = plt.subplots( figsize=figsize )
    
    if col_map is None: col_map='RdBu_r'
    
    if cbar_levels is None: cbar_levels=51
    
    cs = ax.contourf( xx, yy, v, cmap=col_map, levels=cbar_levels, extend='both')
    
    if sonic: 
        with open(sonic,'rb') as f:
            lines = pickle.load( f )
        
        for line in lines:
            x_sonic = line[:,0]
            y_sonic = line[:,1]
            ax.plot(x_sonic,y_sonic,'magenta',linewidth=1.5)


    if separation:
        with open(separation,'rb') as f:
            lines = pickle.load( f )
            
        for line in lines:
            x_sep = line[:,0]
            y_sep = line[:,1]
            ax.plot(x_sep,y_sep,'red',linewidth=1.5)


    if boundary is not None:
        with open( boundary,'rb') as f:
            lines = pickle.load( f )
            
        for line in lines:
            x_bou = line[:,0]
            y_bou = line[:,1]
            ax.plot(x_bou,y_bou,'black',linewidth=1.5)
            
    if shockshape is not None:
        with open(shockshape, 'rb') as f:
            lines = pickle.load( f )
        
        for line in lines:
            x_shock = line[:,0]
            y_shock = line[:,1]
            ax.plot(x_shock,y_shock,'black',linewidth=1.5)
            
    if DS is not None:
        with open(DS, 'rb') as f:
            lines = pickle.load( f )
        
        for line in lines:
            x_DS = line[:,0]
            y_DS = line[:,1]
            ax.plot(x_DS,y_DS,'black',linewidth=0.8)
            
    
    if wall:
        x_wall = np.array([np.min(xx),np.max(xx)])
        y_wall = np.array([-0.1, -0.1])
#        ax.plot(z_wall,y_wall,'black')
        ax.fill_between(x_wall,-0.2,y_wall,color='grey')
        ax.plot(x_wall,np.array([0.0,0.0]),'--',linewidth=0.2)
        
    if tag:
        ax.text( tag_loc[0],tag_loc[1],
                 tag,
                 va='center',
                 fontsize=40,
                 zorder=101)
    
    if x_lim is not None: ax.set_xlim(x_lim)
    if y_lim is not None: ax.set_ylim(y_lim)   
    
    if not pure:
        cbar = plt.colorbar(cs,
                            orientation='horizontal', 
                            shrink=0.5,
                            location='top',
                            aspect=15,
                            ticks=cbar_ticks,
                            pad=0.1)
        
        cbar.outline.set_linewidth(2)
        cbar.ax.set_xlabel(cbar_label,
                           rotation='horizontal',
                           labelpad=20)
        cbar.ax.tick_params( axis='x',
                             direction='in',
                             bottom = True, top=False,
                             labeltop=False,labelbottom=True,
                             length=20.0,
                             width=2.0)
        
        
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(2)
        ax.tick_params(axis='both', length=15, width=1.5, pad=10)
        ax.set_xlabel( r"$(x-x_{imp})/\delta_0$")
        ax.set_ylabel( r"$y/\delta_0$" )

    # set x,y unit length equal
    plt.gca().set_aspect('equal', adjustable='box')

    if filename:
        
        if pure: 
            filename+='_pure'
            fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
            
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
                      tag=None,
                      wall=True, 
                      sonic=True,
                      separation=False,
                      vectors=None, 
                      figsize=None, 
                      filename=None,
                      col_map=None,
                      cbar_label=None,
                      cbar_levels=None,
                      cbar_ticks=None,
                      title=None,
                      arrow=False,
                      extreme_loc=None,
                      x_lim=None,
                      y_lim=None,
                      pure=False):
    
    """
    zz,yy are required to be equally spaced by streamplot
    """

    if figsize is None:
        figsize = (15,9)
    
    fig, ax = plt.subplots( figsize=figsize )

    if col_map is None: col_map='RdBu_r'
    
    if cbar_levels is None: cbar_levels=51

    # extend : fill the region where value exceed min and max
    cs = ax.contourf(zz, yy, v, cmap=col_map, levels=cbar_levels, extend='both')
    
    ax.set_xlim([np.min(zz),np.max(zz)])
    ax.set_ylim([-0.1,1.1])
    
    if pure:
        textfontsize = 80
    else:
        textfontsize = 50

    if sonic: 
        with open('soniclines.pkl','rb') as f:
            lines = pickle.load( f )
        
        for line in lines:
            x_sonic = line[:,0]
            y_sonic = line[:,1]
            ax.plot(x_sonic,y_sonic,'lime',linewidth=5.0)
    
    
    if separation:
        with open('separationlines.pkl','rb') as f:
            lines = pickle.load( f )
            
        for line in lines:
            x_sep = line[:,0]
            y_sep = line[:,1]
            ax.plot(x_sep,y_sep,'red',linewidth=0.8)
    
    if vectors is not None:
        
        if arrow:
            ax.quiver( zz[::8,::6],yy[::8,::6], 
                    vectors[0][::8,::6], vectors[1][::8,::6],
                    width=0.0020,
                    angles='xy',
                    scale_units='xy', 
                    scale = 80,
                    headwidth=5,
                    headlength=7)
            
            rect = matplotlib.patches.Rectangle( (-0.85,0.875), 0.60,0.15, 
                                                facecolor='white',
                                                alpha=0.8 ) # transpanrency
            ax.add_patch(rect)
            
            ax.quiver(-0.84,0.95,
                    15.21,0.0,    # 3% of u_infty
                    width=0.0020,
                    angles='xy',
                    scale_units='xy',
                    scale = 80,
                    headwidth=5,
                    headlength=7)
            
            ax.text( -0.61,0.94,
                     r"$3\%u_{\infty}$",
                     va='center',
                     fontsize=textfontsize)
        else:   # streamline
#            speed = np.sqrt( vectors[0]**2+ vectors[1]**2 )
#            lw = speed/speed.max()
            ax.streamplot( zz,yy, 
                           vectors[0], vectors[1],
                           density=2.0,            # linewidth=lw,
                           color='black',
                           linewidth=2.0)
        
    
    if wall:
        df_wall = pd.read_csv( 'wall_X.dat', delimiter=r'\s+' )
        z_wall = np.array(df_wall['x'])
        y_wall = np.array(df_wall['y'])
#        ax.plot(z_wall,y_wall,'black')
        ax.fill_between(z_wall,-0.1,y_wall,color='grey',zorder=10)

    if tag:
#        rect2 = matplotlib.patches.Rectangle( (0.20,0.85), 0.7,0.20, 
#                                              facecolor='white',
#                                              alpha=0.8,
#                                              zorder=100 ) # transpanrency
#        ax.add_patch( rect2 )
        
        ax.text( 0.1, 0.95,
                 tag,
                 va='center',
                 fontsize=textfontsize,
                 zorder=101,
                 bbox={"fc":"white","alpha":0.8,"ec":"None"})    
    
    # plot a cross at extreme_loc
    if extreme_loc:
        for loc in extreme_loc:
            ax.plot( loc[0],loc[1], 'x', color='black', 
                     markersize=20, markeredgewidth=4.0) 
    
    if not pure: 
           
        if title is not None:
            ax.text(0.5,1.1,title,fontsize=20,transform=ax.transAxes)
        
        # set colorbar
        
        cbar = plt.colorbar(cs,
                            orientation='vertical', 
                            shrink=0.8,
                            location='left',
                            aspect=10,
                            ticks=cbar_ticks,
                            pad=0.15)
        
        cbar.outline.set_linewidth(2)
        
#        cbar.ax.set_xlabel( cbar_label,
#                            fontsize=50,
#                            labelpad=25,
#                            loc='left')
        
        cbar.ax.tick_params( labelsize=50, 
                             direction='in',
                             left=True,right=False,
                             labelleft=True,labelright=False,
                             length=20.0,
                             width=2.0)

        # set ticks on the main xy axises
        
        ax.tick_params(axis='x',labelsize=30,pad=10)
        ax.tick_params(axis='y',labelsize=30,pad=10)

        plt.gca().set_aspect('equal', adjustable='box')

    # set axises stay with contour and x,y unit length equal
    
#    plt.axis('scaled')

    ax.set_xlabel(r"$z/\delta_0$")
    ax.set_ylabel(r"$y_s/\delta_0$")
    ax.tick_params(axis='both', length=10, width=2, pad=10)
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)

    if filename:
        if pure: 
            filename+='_pure'
            fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
            
        plt.savefig( 'figx_'+filename )
        
        print(f"{filename} is output.\n")

    plt.close() 


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/11/08  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_stream_function( zz, yy, v,
                          tag=None,
                          wall=True,
                          figsize=None,
                          filename=None,
                          clevels=None,
                          title=None,
                          extreme_loc=None,
                          fmt='.png',
                          pure=False):
    """
    zz,yy  are required to be equally spaced by streamplot
    """
    
    if figsize is None:
        figsize = (15,9)
    
    fig, ax = plt.subplots( figsize=figsize )
    
    if clevels is None: clevels=21
    
    else:
        clevels_pos = clevels[ clevels>0 ]
        clevels_neg = clevels[ clevels<0 ]
    
    cs_pos = ax.contour( zz, yy, v, 
                         colors='black', 
                         levels=clevels_pos, 
                         linestyles='-',
                         linewidths=4.0)
    cs_neg = ax.contour( zz, yy, v, 
                         colors='black', 
                         levels=clevels_neg, 
                         linestyles='--',
                         linewidths=4.0)

    ax.set_xlim([np.min(zz),np.max(zz)])
    ax.set_ylim([-0.1,1.1])
    
    if pure:
        textfontsize = 80
    else:
        textfontsize = 50
    
    if wall:
        df_wall = pd.read_csv( 'wall_X.dat', delimiter=r'\s+' )
        z_wall = np.array(df_wall['x'])
        y_wall = np.array(df_wall['y'])
        ax.fill_between(z_wall,-0.1,y_wall,color='grey',zorder=10)
    
    if tag:
        ax.text( 0.1, 0.95,
                 tag,
                 va='center',
                 fontsize=textfontsize,
                 zorder=101,
                 bbox={"fc":"white","alpha":0.8,"ec":"None"})
    
    # plot a cross at extreme_loc
    if extreme_loc:
        for loc in extreme_loc:
            ax.plot( loc[0],loc[1], 'x', color='black', 
                     markersize=20, markeredgewidth=4.0)
    
    if not pure:
        
        if title is not None:
            ax.text(0.5,1.1,title,fontsize=20,transform=ax.transAxes)
        
        # set ticks on the main xy axises
        
        ax.tick_params(axis='x',labelsize=30,pad=10)
        ax.tick_params(axis='y',labelsize=30,pad=10)
        
        plt.gca().set_aspect('equal', adjustable='box')
        
    
    if filename:
        if pure:
            filename+='_pure'
            fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        
        filename += fmt
        plt.savefig( 'figx_'+filename )
        
        print(f"{filename} is output.\n")
    
    plt.close()



# ----------------------------------------------------------------------
# >>> plot friction projection                                    (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/26  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_wall_projection( xx, zz, v,
                          separation=None,
                          figsize=None,
                          filename=None,
                          col_map=None,
                          cbar_label=None,
                          cbar_levels=None,
                          cbar_ticks=None,
                          label=None,
                          title=None,
                          extend='both',
                          use_norm=False,
                          pure=False,):
    
    """
    only applicable to the ridge type smooth wall case setting.
    Range in z is manually set [-2.0,2.0]
    """

    if figsize is None:
        if pure: figsize=(30,4) 
        else:    figsize = (25,10)

    fig, ax = plt.subplots( figsize=figsize )

    if col_map is None: col_map='coolwarm'

    if cbar_levels is None: cbar_levels=51

    zz[0,:] =-2.0   # for the convenience of plotting
    zz[-1,:]= 2.0
    
    cs = ax.contourf( xx, zz, v,
                      levels=cbar_levels,
                      cmap=col_map,
                      extend=extend,
                      norm=colors.CenteredNorm() if use_norm else None)
    
    if separation is not None:
        
        with open(separation,'rb') as f:
            
            lines = pickle.load( f )
            for line in lines:
                x_sep = line[:,0]
                z_sep = line[:,1]
                ax.plot(x_sep,z_sep,'black',linewidth=1.2)
                
    ax.set_ylim([-2.0,2.0])
    ax.set_xlim([-20.0,10.0])
    
    if not pure:
        cbar = plt.colorbar(cs,orientation='horizontal', shrink=0.7, ticks=cbar_ticks)
        cbar.ax.set_ylabel(cbar_label,rotation=0,fontsize=30,labelpad=150,loc='bottom')
        cbar.ax.tick_params(labelsize=30)
        cbar.ax.set_position([0.4,0.7,0.3,0.1])
        # cbar.ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        
        ax.minorticks_on()
        ax.tick_params( which='major',
                        axis='both',
                        direction='out',
                        length=15,
                        width=2)
        ax.tick_params( which='minor',
                        axis='both', 
                        direction='out',
                        length=10,
                        width=1)

        ax.tick_params(axis='x',labelsize=30,pad=10)
        ax.tick_params(axis='y',labelsize=30,pad=10)
        
        ax.set_xlabel(r"$(x-x_{imp})/\delta_0$")
        ax.set_ylabel(r"$z/\delta_0$")
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(1.5)
        
        ax.text(5,-1,label)
        
        plt.title(title,fontsize=30,pad=200)
        
        # set axises stay with contour and x,y unit length equal
        
        plt.axis('tight')
        plt.gca().set_aspect('equal', adjustable='box')
    
        

    if filename:
        
        if pure:
            filename +='_pure'
            fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
            
        plt.savefig( filename )
        print(f"{filename} is output.\n")

    plt.close() 


# ----------------------------------------------------------------------
# >>> plot spanwise distribution of wall-projected variables
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/06/05  - created
#
# Desc
#
# ----------------------------------------------------------------------

def plot_spanwise_variables( z, v, ylabel, figname ):
    
    fig, ax = plt.subplots( figsize=(15,10) )
    
    ax.plot( z, v, color='black', marker='o', markersize=10)
    
    ax.set_xlabel(r"$z/\delta_0$")
    ax.set_ylabel(ylabel)
    
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
    
    ax.grid(which='both',visible=True)
    ax.set_xlim([0,2])
#    ax.set_ylim([0,10])
    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(2)
    
    plt.savefig(figname)
    print(f"{figname} is output.\n")
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