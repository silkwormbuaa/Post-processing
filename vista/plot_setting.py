# -*- coding: utf-8 -*-
'''
@File    :   plot_setting.py
@Time    :   2024/09/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   a module defines some common settings for plotting 
             related to pyvista and matplotlib
'''

import pyvista           as     pv
import matplotlib.pyplot as     plt
import matplotlib.ticker as     ticker


# ----------------------------------------------------------------------
# >>> cpos_callback                                              (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/09/02  - created
#
# Desc 
#   in the pyvista.RenderWindowInteractor, add a callback function to
# ----------------------------------------------------------------------

def cpos_callback( p:pv.Plotter ):
    
    """
    In the interactive window, press 'p' to print the camera position in the console
    
    Des: in the pyvista.RenderWindowInteractor, add a callback function to
    """
    
    def cpos_print( p:pv.Plotter ):
        
        pfv = p.camera_position.to_list()
        pos = [f"{float(num):5.2f}" for num in pfv[0]]
        foc = [f"{float(num):5.2f}" for num in pfv[1]]
        vup = [f"{float(num):5.2f}" for num in pfv[2]]
        print(f"Camera position: {pos}")
        print(f"Camera focus: {foc}")
        print(f"Camera viewup: {vup}")
        print("======")
        
        return
    
    p.add_key_event( "p", lambda: cpos_print(p) )

    return


# ----------------------------------------------------------------------
# >>> set plt rcparams                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/09/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def set_plt_rcparams(preamble='stix'):

    plt.rcParams["text.usetex"] = True
    
    if preamble == 'stix':
        plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
    elif preamble == 'amssymb':
        plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
        
    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size']   = 40


# ----------------------------------------------------------------------
# >>> set plt style                                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/09/09  - created
#
# Desc
#
# ----------------------------------------------------------------------

def set_plt_style( case=None, ax=None, fig=None ):    
    # pressure fluctuation in streamwise direction

# =============================================================================

    if case == 'pressure_fluctuation':
        
        # set figure size
        fig.set_size_inches(15,8)

        # set contrained layout
        fig.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.2)
        
        # set the ranges of x,y axis
        ax.set_xlim([-15,10])
        ax.set_ylim([0.01,0.09])

        # set label and distance
        ax.set_xlabel(r"$(x-x_{imp})/\delta_0$", labelpad=-5 )  
        ax.tick_params(axis='x', pad=15)

        ax.set_ylabel(r"$\sqrt{\langle p'p' \rangle}/p_{\infty}$" )
        ax.tick_params(axis='y', pad=10)

        # set ticks
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))

        ax.minorticks_on()
        ax.tick_params( which='major',
                        axis='both',
                        direction='in',
                        length=15,
                        width=2 )
        ax.tick_params( which='minor',
                        axis='both', 
                        direction='in',
                        length=10,
                        width=1 )

        # set the bounding box of axes
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(3)

# =============================================================================

    elif case == 'wall_pressure':
        
        # set figure size
        fig.set_size_inches(15,8)

        # set contrained layout
        fig.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.2)
  
        # set the range of x,y axis 
        ax.set_xlim([-15,10])
        ax.set_ylim([0.8, 2.5])
        
        # set x,y label and distance
        ax.set_xlabel("$(x-x_{imp})/\delta_0$", labelpad=-5 )
        ax.tick_params(axis='x', pad=15)
        
        ax.set_ylabel("$<p_w>/p_{\infty}$")
        ax.tick_params(axis='y', pad=10)

        # set ticks
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))

        ax.minorticks_on()
        ax.tick_params( which='major', 
                        axis='both',
                        direction='in',
                        length=15,
                        width=2 )
        ax.tick_params( which='minor',
                        axis='both', 
                        direction='in',
                        length=10,
                        width=1 )
        
        # set the bounding box of axes
        ax.spines[:].set_color('black')
        ax.spines[:].set_linewidth(3)  
        
        
# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/09/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    pass



# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/09/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()