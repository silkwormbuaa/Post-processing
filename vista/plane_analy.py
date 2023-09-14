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

import pickle

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
    
    lines = []
    
    fig, ax = plt.subplots(figsize=(10, 4))

    cs = ax.contour( xx, yy, mach, levels=[1.0] )

    for isoline in cs.collections[0].get_paths():
        line = isoline.vertices
        lines.append( line )
        
    plt.close()
    
    with open('soniclines.pkl','wb') as f:
        
        pickle.dump( lines, f )
        


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
    
    return: pickle isolines list into 'separationlines.pkl'
    """
    
    lines = []
    
    fig, ax = plt.subplots(figsize=(20, 8))
    
    cs = ax.contour(xx, yy, u, levels=[0.0] )
    
    for isoline in cs.collections[0].get_paths():
        line = isoline.vertices
        lines.append( line )
        
    plt.close()
    
    with open('separationlines.pkl','wb') as f:
        
        pickle.dump( lines, f)



# ----------------------------------------------------------------------
# >>> Save isolines                                               (Nr.)
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

def save_isolines( xx, yy, v, value:float, file ):

    """
    xx,yy: 2d numpy array storing coordinates
    v : 2d numpy array of target variable
    
    return: pickle isolines list into file
    """
    
    lines = []
    
    fig, ax = plt.subplots(figsize=(20, 8))
    
    cs = ax.contour(xx, yy, v, levels=[0.0] )
    
    for isoline in cs.collections[0].get_paths():
        line = isoline.vertices
        lines.append( line )
        
    plt.close()
    
    with open( file,'wb') as f:
        
        pickle.dump( lines, f)
        


# ----------------------------------------------------------------------
# >>> shift coordinates                                          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/09/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

def shift_coordinates( df:pd.DataFrame, delta, H, H_md, x_imp ):
    
    """
    df    : dataframe of df_slice
    delta : boundary thickness
    H     : height of ridge
    H_md  : melt-down height
    x_imp : x coordinate of impingement point
    
    return : df with xs, ys, zs (if x,y,z exist)
    """
    
    vars = df.columns
    
    if 'x' in vars: 
        xs = ( np.array( df['x'] ) - x_imp ) / delta
        df['xs'] = xs
    
    if 'y' in vars: 
        ys = ( np.array( df['y'] ) + H - H_md ) / delta
        y_scale = np.array( df['y'] ) / delta
        df['ys'] = ys
        df['y_scale'] = y_scale
        
    if 'z' in vars:
        zs = np.array( df['z'] ) / delta
        df['zs'] = zs
    
    return df



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
