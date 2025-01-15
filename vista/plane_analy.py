# -*- coding: utf-8 -*-
'''
@File    :   2d_analysis.py
@Time    :   2023/09/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   2D data analysis 
'''

import os
import pickle
import numpy             as     np
import pandas            as     pd
import pyvista           as     pv
import matplotlib.pyplot as     plt

from   .snapshot         import Snapshot
from   .grid             import GridData

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

def save_sonic_line( xx, yy, mach, out_file='soniclines.pkl' ):
    
    """
    xx,yy : 2d numpy array storing coordinates
    mach  : 2d numpy array of mach number field
    out_file : output file name
    
    content of out_file: 
    1. list of sonic lines(each line is a 2D array of coordinates )
    2. list of mean height of sonic lines
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
    
    with open(out_file,'wb') as f:
        
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

def save_isolines( xx, yy, v, value:float, file, clip=False ):

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
        
        if clip and all( line[:,1] < 0.5 ):
            continue
        
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
    array  = np.array( array ).reshape( N, int(n0/N), n1 )
    array  = np.mean( array, axis=0 )
    
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

def compute_DS( grad_rho, min=None, max=None ):
    
    """
    grad_rho : array of density gradient
    min, max : min and max value of grad_rho, if not given, will be computed from grad_rho
    
    return : same shape array of DS. 
    Refer to Wu and Martin(2007) and Guo(2022) for more details.
    """
    
    DS = np.zeros_like( grad_rho )
    
    if min is None:  min = np.min( grad_rho )
    if max is None:  max = np.max( grad_rho )
    
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
    integrate stream function (should only be applied to 2D flows or 3D axisymmetric flows)
    
    rho: 2D array of density, shape(ny,nz)
    rho_w: reference density at wall
    w_favre: 2D array of favre averaged horizontal velocity
    v_favre: 2D array of favre averaged vertical velocity
    
    zz: z coordinate
    yy: y coordinate
    
    """
    
    psi = np.zeros_like( rho )
    
# - integrate w in -y direction
    
    dy = -np.diff(y)
    
    for j in range(len(y)-1,0,-1):
        psi[j-1,-1] = psi[j,-1] + (dy[j-1] * 0.5*(w_favre[j,-1]*rho[j,-1]+w_favre[j-1,-1]*rho[j-1,-1]))

# - integrate -v in -z direction
    
    dz = -np.diff(z)
    
    for i in range(len(z)-1,0,-1):
        psi[:,i-1] = psi[:,i] - (dz[i-1] * 0.5*(v_favre[:,i]*rho[:,i]+v_favre[:,i-1]*rho[:,i-1])) 
    
# - normalize

    psi = psi / rho_w
    
    return psi 



# ----------------------------------------------------------------------
# >>> Interpolate data into a 2D plane                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/01/12  - created
#
# Desc
#
# ----------------------------------------------------------------------

def pv_interpolate( source:pv.MultiBlock, vars, mesh_vectors ) -> pd.DataFrame:
    
    """
    source: data source, a pyvista MultiBlock object.
    vars:   list of variables
    mesh_vectors: [px,py,pz], 3D, the value of cross-plane direction is 0.
    
    Return: a data frame. Data should be reshaped as (npj,npi).
     
    interpolate data into a uniform or rectilinear cartesian grid.
    """
    
    px = np.array( mesh_vectors[0] ); npx = len(px)
    py = np.array( mesh_vectors[1] ); npy = len(py)
    pz = np.array( mesh_vectors[2] ); npz = len(pz)
    
    grid = pv.RectilinearGrid( px, py, pz )
    
    interped = grid.sample( source )

    # using 'ij' indexing will return matrix with size (npx,npy,npz)
    
    x,y,z = np.meshgrid( px, py, pz, indexing='ij' )
    
    # since data should be reshaped as (npj,npi) if you want to plot streamline
    # with matplotlib, so we should adjust the order of coordinates, which 
    # originally in the order of (npx,npy,npz) (C-format).
    
    if   npx == 1: pass
    elif npy == 1:
        z = z.reshape( npx, npz ).T
        x = x.reshape( npx, npz ).T
    elif npz == 1:
        y = y.reshape( npx, npy ).T
        x = x.reshape( npx, npy ).T
    else: 
        raise ValueError("Only support 2D data interpolation!")
    
    df = pd.DataFrame( {'x':x.ravel(),'y':y.ravel(),'z':z.ravel()} )
    
    # pv.RectlinearData.point_data are ordered as ( npz*npy*npx, n_vars)
    
    if   npx == 1:  
        for var in vars:
            df[var] = np.array(interped.point_data[var]).reshape(npz,npy).T.ravel()
    
    elif npy == 1 or npz ==1:  
        for var in vars:
            df[var] = np.array(interped.point_data[var])
    
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

    snapfile = '/home/wencan/temp/220927/snapshots/snapshot_01327116/snapshot_Z_001.bin' 
    gridfile = '/home/wencan/temp/220927/results/inca_grid.bin'
    
    os.chdir('/home/wencan/temp/test')
    
    snap = Snapshot( snapfile )
    snap.read_snapshot()
    
    grid = GridData( gridfile )
    grid.read_grid()
    
    snap.grid3d = grid
    dataset = pv.MultiBlock( snap.create_vtk_multiblock() )
    
    print( dataset )

    px = np.linspace(0,   20, 101, endpoint=True)
    py = np.linspace(0.0,  10.4, 201, endpoint=True)
    pz = np.array( [0.0] )
    
    df = pv_interpolate( dataset, ['u','v','w','T'], [px,py,pz] )
    
    print( df )
    
    u = np.array( df['u'] ).reshape( 201, 101 )
    v = np.array( df['v'] ).reshape( 201, 101 )
    w = np.array( df['w'] ).reshape( 201, 101 )
    z = np.array( df['z'] ).reshape( 201, 101 )
    x = np.array( df['x'] ).reshape( 201, 101 )  
    y = np.array( df['y'] ).reshape( 201, 101 )

    
    # print( len(interped.point_data["u"]) )

    fig, ax = plt.subplots( figsize = [12,8] )
    
    cs = ax.contourf( x,y, u)
    
    
    ax.streamplot( x,y, 
                   u,v,
                   density=2,
                   color='black')
    
    plt.savefig("test.png")

    plt.close()    

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
