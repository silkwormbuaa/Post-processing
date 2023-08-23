# -*- coding: utf-8 -*-
'''
@File    :   geometry.py
@Time    :   2023/08/16 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import sys

import os

source_dir = os.path.realpath(__file__).split('vista')[0] 
sys.path.append( source_dir + 'vista/lib/form' )

print(source_dir + 'lib/exec/')

import numpy             as     np

from   stl               import mesh

import time

from   .math_opr         import unitize_L2

from   .lib.form         import geo



# ----------------------------------------------------------------------
# >>> Read STL file                                              (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/16  - created
#
# Desc
#
#   - read a stl file using numpy-stl library
#
# ----------------------------------------------------------------------

def read_stl( filename ):
    
    msh = mesh.Mesh.from_file( filename )
    
    # unitize the normals (numpy-stl returns normals directly by calculation)
    
    for i, norm in enumerate(msh.normals):
    
        msh.normals[i] = unitize_L2( norm )
    
    # calculate the number of triangles
    
    msh.n_tri = len( msh.normals )
        
    # store the bounding box coordinates
    
    msh.minx = min( np.min(msh.v0[:,0]),
                    np.min(msh.v1[:,0]),
                    np.min(msh.v2[:,0]) )

    msh.maxx = max( np.max(msh.v0[:,0]),
                    np.max(msh.v1[:,0]),
                    np.max(msh.v2[:,0]) )
    
    msh.miny = min( np.min(msh.v0[:,1]),
                    np.min(msh.v1[:,1]),
                    np.min(msh.v2[:,1]) )

    msh.maxy = max( np.max(msh.v0[:,1]),
                    np.max(msh.v1[:,1]),
                    np.max(msh.v2[:,1]) )
    
    msh.minz = min( np.min(msh.v0[:,2]),
                    np.min(msh.v1[:,2]),
                    np.min(msh.v2[:,2]) )

    msh.maxz = max( np.max(msh.v0[:,2]),
                    np.max(msh.v1[:,2]),
                    np.max(msh.v2[:,2]) )
    
    return msh



# ----------------------------------------------------------------------
# >>> Point in triangle b                                         (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/16  - created
#
# Desc
#   - Adapted from inca.inca_stl_geo.POINT_IN_TRIANGLE_B
# ----------------------------------------------------------------------

def point_in_triangle_b( p, p1, p2, p3 ):
    
#-- initialize status for quick return
    
    flag = False

#-- take p1 as origin

    v0 = p3 - p1
    v1 = p2 - p1
    v2 = p  - p1
    
#-- compute dot products
    
    dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2] # v0 o v0
    dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2] # v0 o v1
    dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2] # v0 o v2
    dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] # v1 o v1
    dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] # v1 o v2
    
#-- compute barycentric coordinates

    d = dot00 * dot11 - dot01 * dot01
    u = dot11 * dot02 - dot01 * dot12
    v = dot00 * dot12 - dot01 * dot02
    
#-- evaluate

    if( d >= 0.0 ):
        flag = ( u >= 0.0 and v >= 0.0 and u + v <= d )
    elif( d < 0.0 ):
        flag = ( u <= 0.0 and v <= 0.0 and u + v >= d )
    
    
    return flag



# ----------------------------------------------------------------------
# >>> ray_tracer                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/16  - created
#
# Desc
#
#   - using ray-tracing method to assign inside/outside cells
#   - x_init and signum_init should correspond. 
#     Recommend x_init to be outside ib and signum_init to be 1.0.
#   - return value is vol_fra[i,j,k] with inside 0.0 and outside 1.0.
# ----------------------------------------------------------------------

def ray_tracer( ibmsh, g, x_init, signum_init, buff = 3, verbose=False ):

#-- parameters

    C13 = np.float64(1/3)
    EPS = sys.float_info.epsilon * 2.0
    
#-- Initialization
    
    # volume fraction array
    
    vol_fra = np.zeros( shape=(g.nx+buff*2,g.ny+buff*2,g.nz+buff*2), dtype='f' )
    
    # signum_init represent the signum at x_init
    # signum must be a numpy array to pass through Fortran-python interface
    
    signum = np.array(signum_init,dtype=np.float64)
    
    # reference point
    
    x_ref = x_init    

    # number of grid point in each direction
    
    npx = g.nx + buff*2
    npy = g.ny + buff*2
    npz = g.nz + buff*2
    
    # lam_list
    
    lam_list = np.zeros(1000)


#-- Move reference location as closely as possible to IB

    if ( (x_ref[0] < ibmsh.minx) or (x_ref[0] > ibmsh.maxx) or
         (x_ref[1] < ibmsh.miny) or (x_ref[1] > ibmsh.maxy) or
         (x_ref[2] < ibmsh.minz) or (x_ref[2] > ibmsh.maxz) ):
        
        # define ray from x_ref to block center
        
        x  = np.array( [ g.gx[npx//2], g.gy[npy//2], g.gz[npz//2] ] )
        dx = x - np.array(x_ref)
        
        # calculate intersections
        
        lam_box = [ ( ibmsh.minx - x_ref[0] ) / dx[0],
                    ( ibmsh.maxx - x_ref[0] ) / dx[0],
                    ( ibmsh.miny - x_ref[1] ) / dx[1],
                    ( ibmsh.maxy - x_ref[1] ) / dx[1],
                    ( ibmsh.minz - x_ref[2] ) / dx[2],
                    ( ibmsh.maxz - x_ref[2] ) / dx[2] ]
        
        lam_min = max( min( lam_box[0], lam_box[1] ),
                       min( lam_box[2], lam_box[3] ),
                       min( lam_box[4], lam_box[5] ) )
        
        lam_max = min( max( lam_box[0], lam_box[1] ),
                       max( lam_box[2], lam_box[3] ),
                       max( lam_box[4], lam_box[5] ) )
        
        # cast ray( x_ref -> x ) on IB bounding box
        
        if ( lam_min > lam_max ):
            
            # line does not cut bounding box 
            # set new reference point to block center
            x_ref = x_ref + 0.999 * dx
            
        elif( (lam_max < 0.0) or (lam_min > 1.0) ):
            
            # line does cut bounding box, but outside ray interval 
            # set new reference point to block center
            x_ref = x_ref + 0.999*dx
            
        elif( lam_min < 0.0 ):
            
            # ray exits bounding box at lam_max 
            # (this should not happen unless ray starts in bounding box)
            # keep current x_ref
            pass
        
        else:
            
            # ray enters bounding box at lam_min
            # set new reference point to lam_min (where ray enters box)
            x_ref = x_ref + 0.999 * lam_min * dx
    

#-- Loop over all cells in the block

    for k in range( npz ):
        for j in range( npy ):
            for i in range( npx ):

                # define the ray from reference point to grid point
                x  = np.array( [ g.gx[i], g.gy[j], g.gz[k] ] )
                dx = x - x_ref
            
                # overlap of ray ('s bounding box) with ib bounding box
                
                x_max = np.maximum( x_ref, x )
                x_min = np.minimum( x_ref, x )
                cover = [ min(ibmsh.maxx,x_max[0]) - max(ibmsh.minx,x_min[0]),
                          min(ibmsh.maxx,x_max[0]) - max(ibmsh.minx,x_min[0]),
                          min(ibmsh.maxx,x_max[0]) - max(ibmsh.minx,x_min[0]) ]
                
                # all elements in cover >= 0, ray is possible intersecting ib
                
                if ( all( item >= 0.0 for item in cover ) ):             

                    # loop over all ibmsh triangles

                    cnt = 0
                    
                    for i_tri in range( ibmsh.n_tri ):
                        
                        p0 = ibmsh.v0[i_tri]
                        p1 = ibmsh.v1[i_tri]
                        p2 = ibmsh.v2[i_tri]
                        normal = ibmsh.normals[i_tri]

                        # compute overlap of triangle and block 
                        
                        corner_up = [max(u,v,w) for u,v,w in zip(p0,p1,p2)]
                        corner_dn = [min(u,v,w) for u,v,w in zip(p0,p1,p2)]

                        cover = np.minimum( np.array(corner_up), x_max ) \
                              - np.maximum( np.array(corner_dn), x_min )
                        
                        if ( all( item >= 0.0 for item in cover) ):
                            
                            # ray's projection length on normal direction
                            
                            v_n = np.dot( normal, dx )
                            
                            if ( v_n != 0 ): # ray not parallel to triangle
                                
                                tri_cen = C13 * ( p0 + p1 + p2 )
                                
                                lam = np.dot( normal, (tri_cen-x_ref) ) / v_n
                                
                                if ( 0.0 <= lam <= 1.0 ):
                                    
                                    # check if this intersection point 
                                    # was found before
                                    
                                    tol = sum(np.maximum(abs(tri_cen),
                                                         abs(x_ref))) * EPS
                                    
                                    new = True
                                    
                                    for ii in range( cnt ):
                                        if ( abs(lam - lam_list[ii]) < tol ):
                                            new = False ; break
                                        
                                    # if triangle intersects ray in a new point
                                    if new:
                                        
                                        ip = x_ref + dx * lam                                        
                                        
                                        if point_in_triangle_b(ip,p0,p1,p2):
                                            
                                            # found valid intersection
                                            
                                            lam_list[cnt] = lam
                                            cnt = cnt + 1
                                            if cnt == 1000: 
                                                cnt = 998     # buffer overflow    
                                    # new
                                # intersection within ray range
                            # triangle not parallel to ray
                        # triangle (or partly) within block
                    # loop over all triangles
                    
                    if cnt%2 == 1 : 
                        
                        # odd intersection points, x and x_ref have opposite sig
                        # here just distinguish inside/outside ib, later, 
                        # interface cells vol_fra will be overwritten.
                        signum = -signum
            
#                    geo.ray_tracer_core(x,x_ref,ibmsh.normals,
#                                        ibmsh.v0,ibmsh.v1,ibmsh.v2,
#                                        signum)
                else:
                    
                    # ray and ib won't intersect, keep signum
                    pass
                
                # give signum to vol_frc at [i,j,k], and take [i,j,k] as 
                # reference point for next cell
                
                vol_fra[i,j,k] = signum
                x_ref = x      
                      
            # end i
        # end j
    # end k              
            
    vol_fra[vol_fra == -1.0] = 0.0
    
    return vol_fra



# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    filename = '/home/wencanwu/testdir/sphere.stl'
    
    mesh = read_stl( filename )
    
    print( mesh.normals[0] )
    
    print( mesh.v0[0],mesh.v1[0],mesh.v2[0] )
    
    print(mesh.v0.shape)
    
    print( mesh.minx )
    print( mesh.maxx )
    print( mesh.miny )
    print( mesh.maxy )
    print( mesh.minz )
    print( mesh.maxz )
    
    print(f"number of triangles: {mesh.n_tri}")

    p = np.array([0.0,0.5,0.0])
    p1 = np.array([1.0,0.0,0.0])
    p2 = np.array([-1.0,0.0,0.0])
    p3 = np.array([0.0,2.0,0.0])
    
    print( geo.point_in_triangle.__doc__)
    
    flag = bool( geo.point_in_triangle(p,p1,p2,p3) )
    
    print( flag )


# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/08/16  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()