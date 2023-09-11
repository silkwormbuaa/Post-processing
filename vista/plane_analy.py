# -*- coding: utf-8 -*-
'''
@File    :   2d_analysis.py
@Time    :   2023/09/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   2D data analysis 
'''

import numpy             as     np

import pandas            as     pd

import matplotlib.pyplot as     plt

# ----------------------------------------------------------------------
# >>> Save sonic line                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

def save_sonic_line( xx, yy, mach ):
    
    """
    xx,yy: 2d numpy array storing coordinates
    mach : 2d numpy array 
    """
    
    header = 'x      y'
    xycor = np.empty(shape=[0, 2])
    fig, ax = plt.subplots(figsize=(10, 4))
    cs = ax.contour(xx, yy, mach, levels=[1.0], linewidths=1.5, colors="k")
    for isoline in cs.collections[0].get_paths():
        xy = isoline.vertices
        xycor = np.vstack((xycor, xy))
        
    plt.close()
        
    np.savetxt(
        "SeprationLine.dat",
        xycor,
        fmt="%15.7f",
        delimiter=" ",
        comments="",
        header=header,
    )
    


# ----------------------------------------------------------------------
# >>> Save separation line                                       (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

def save_separation_line( xx, yy, u ):

    """
    xx,yy: 2d numpy array storing coordinates
    u : 2d numpy array of streamwise velocity
    """
    
    header = 'x      y'
    xycor = np.empty(shape=[0, 2])
    fig, ax = plt.subplots(figsize=(10, 4))
    cs = ax.contour(xx, yy, u, levels=[1.0], linewidths=1.5, colors="k")
    for isoline in cs.collections[0].get_paths():
        xy = isoline.vertices
        xycor = np.vstack((xycor, xy))
    
    plt.close()
    
    np.savetxt(
        "SonicLine.dat",
        xycor,
        fmt="%15.7f",
        delimiter=" ",
        comments="",
        header=header,
    )

# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/11  - created
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
# 2023/09/11  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
