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
    lines: list of sonic lines(each line is a 2D array of coordinates )
    """
    
    lines = []
    mean_heights = []
    
    fig, ax = plt.subplots(figsize=(10, 4))

    cs = ax.contour( xx, yy, mach, levels=[1.0] )

    for isoline in cs.collections[0].get_paths():
        line = isoline.vertices
        lines.append( line )
        
    plt.close()
    
    # compute the average height of sonic line
    
    for i, line in enumerate(lines):
        
        mean_height = np.array( line[:,1] ).mean()
        print(f"Mean height of sonic line [{i}] is {mean_height}.\n")
        mean_heights.append(mean_height)
    
    with open('soniclines.pkl','wb') as f:
        
        pickle.dump( lines, f )
        pickle.dump( mean_heights, f )
        

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
    
    cs = ax.contour(xx, yy, v, levels=[value] )
    
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
# >>> Periodic Average                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

def periodic_average( array, N:int, axis=0, sym=False, antisym=False ):
    
    """
    array    : target array
    N        : number of periods
    axis     : which axis to do periodic average
    
    return: same shape array after periodic average    
    """

    if axis==1:
        array = array.T
        
    n0, n1 = np.array( array ).shape
    array = np.array( array ).reshape( N, int(n0/N), n1 )
    array = np.mean( array, axis=0 )
    
    # symmetrize or anti-symmetrize the array
    
    if sym:
        for i in range( int(n0/N/2) ):
            sym_mean = 0.5*( array[i,:] + array[-i-1,:] )
            array[i,:] = sym_mean
            array[-i-1,:] = sym_mean
    if antisym:
        for i in range( int(n0/N/2) ):
            sym_mean = 0.5*( array[i,:] - array[-i-1,:] )
            array[i,:] = sym_mean
            array[-i-1,:] = -sym_mean
    
    array = np.array( [array]*N ).reshape( n0, n1 )
    
    if axis==1:
        array = array.T
    
    return array


# ----------------------------------------------------------------------
# >>> Compute the separation ratio                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/02  - created
#
# Desc
#
# ----------------------------------------------------------------------

def compute_separation_ratio( array, write=True ):
    
    """
    array: equally spaced 2D array of friction or Cf
    
    return: the ratio of separation area over the whole covered area\n
    write out separation length ratio along z (same size with z)
    """
    
    total_size = array.size
    
    sep_size = np.sum( array < 0.0 )
    
    area_ratio = sep_size/total_size
    
    npz = array.shape[0]
    npx = array.shape[1]
    
    len_ratio = [ np.sum(array[i,:]<0.0)/npx for i in range(npz) ]
    
    if write:
        with open('separation_ratio_area.dat','w') as f:
            f.write( f"{area_ratio:10.5f}" )
        with open('separation_ratio_len.pkl','wb') as f:
            pickle.dump( len_ratio, f )
    
    return area_ratio
    

# ----------------------------------------------------------------------
# >>> compute DS from grad_rho                                     (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def compute_DS( grad_rho ):
    
    """
    grad_rho : array of density gradient
    
    return : same shape array of DS. 
    Refer to Wu and Martin(2007) and Guo(2022) for more details.
    """
    
    DS = np.zeros_like( grad_rho )
    
    min = np.min( grad_rho )
    max = np.max( grad_rho )
    
    DS = 0.8*np.exp( -10.0*(grad_rho-min) / (max-min) )
    
    return DS


# ----------------------------------------------------------------------
# >>> Compute compressible stream function                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/27  - created
#
# Desc
#
# ----------------------------------------------------------------------

def compute_stream_function( rho, rho_w, w_favre, v_favre, z, y ):
    
    """
    integrate stream function 
    
    rho: 2D array of density, shape(ny,nz)
    rho_w: reference density at wall
    w_favre: 2D array of favre averaged horizontal velocity
    v_favre: 2D array of favre averaged vertical velocity
    
    zz: z coordinate
    yy: y coordinate
    
    """
    
    psi = np.zeros_like( rho )
    
    # integrate -v in z direction
    
    dz = np.diff(z)
    
    psi[:,0] = 0.0
    
    for i in range(1,len(z)):
        psi[:,i] = psi[:,i-1] - (dz[i-1] * 0.5*(v_favre[:,i]+v_favre[:,i-1]) 
                                 * rho[:,i])
    
    # integrate w in y direction
    
    dy = np.diff(y)
    
    for i in range(1,len(y)):
        psi[i,:] = psi[i-1,:] + (dy[i-1] * 0.5*(w_favre[i,:]+w_favre[i-1,:])
                                 * rho[i,:])
        
    # normalize
    
    psi = psi / rho_w
    
    return psi 


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
