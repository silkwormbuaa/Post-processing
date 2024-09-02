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

import pyvista           as pv
import matplotlib.pyplot as plt


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
    
    p.add_key_event( "p", cpos_print(p) )

    return

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