#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_zslice.py
@Time    :   2025/04/24 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import pickle
import numpy             as     np
import pyvista           as     pv
import matplotlib.pyplot as     plt
from   scipy.interpolate import Rbf

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.params      import Params
from   vista.directories import Directories
from   vista.tools       import crop_to_rect_map
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams( fontsize=20 )

def main():
    case_dir   = '/home/wencan/temp/220927/' 
    clipbox    = [-15, 10, 0, 10.0, -1, 1]
    cbar_ticks = np.linspace(0.0,2.0,5, endpoint=True)
    
    dataset, x_pfmax = data_preparation( case_dir )
    dataset = dataset.clip_box( clipbox, invert=False )
    
    dirs    = Directories( case_dir )
    os.chdir( dirs.pp_z_average )
    
    pv_visualize( dataset, 'mach', clipbox, cbar_ticks, x_pfmax )


def dividing_streamline( dataset:pv.UnstructuredGrid, point_start:np.ndarray, bbox, step=0.1, forward=True):

    """
    compute the dividing streamline
    """

    u,v = fetch_uv(dataset, point_start)
    pts = [point_start]
    us  = [u]
    vs  = [v]
    
    i = 0
    
    while True:
        
        # - first, try to step forward see if the next point is valid
        
        step_t = step
        point_start_t = step_streamline( pts[-1], us[-1], vs[-1], step=step_t, forward=forward )
        u, v = fetch_uv(dataset, point_start_t)
        
        # - if the next point is invalid, try a larger step size.
        # - the step size is doubled until the next point is valid
        # - or the step size is larger than 1.0.
        
        while u == 0.0 and v == 0.0:
            step_t = step_t*2
            point_start_t = step_streamline( pts[-1], us[-1], vs[-1], step=step_t, forward=forward )
            u,v = fetch_uv(dataset, point_start_t)
            print(f"invalid point in step {i}, point_start_t: {point_start_t}, u: {u}, v: {v}, step_t: {step_t}")
            if step_t >1:
                break
        
        # - if the next point is valid, step forward with the updated step size
        
        point_start = step_streamline( pts[-1], us[-1], vs[-1], step=step_t, forward=forward )
        u, v = fetch_uv(dataset, point_start)
        
        pts.append( point_start )
        us.append( u )
        vs.append( v )
        
        i += 1
        
        if i > 15000:
            break
        
        if point_start[0] < bbox[0] or point_start[0] > bbox[1] or \
           point_start[1] < bbox[2] or point_start[1] > bbox[3]:
            break
        
    return np.array(pts)


def step_streamline( point_start, u, v, step=0.1, forward=True):
    """
    step the streamline
    """

    if not forward:
        step = -step
        
    v_mag = np.sqrt(u**2 + v**2)
    
    x = point_start[0] + step*u/v_mag
    y = point_start[1] + step*v/v_mag
    z = point_start[2]

    return np.array([x,y,z])
    

def fetch_uv( dataset:pv.UnstructuredGrid, point:np.ndarray ):
    """
    fetch the u, v, a at the point
    """
    pt=pv.PolyData(point)
    results = pt.sample( dataset )
    
    return results['u'][0], results['v'][0]
    

def data_preparation( case_dir ):

    dirs     = Directories( case_dir )
    params   = Params( dirs.case_para_file )
    x_pfmax  = params.x_pfmax
    
    # - check if spanwise-averaged vtmb file exists

    filename = dirs.pp_z_average + '/z_average.vtmb'
    
    if not os.path.exists( filename ):
        print(f"Directory {filename} does not exist.")
        print("Please run 'scripts/statistics/spanwise_average.py' first!")
        sys.exit()
    
    else:
        # - read in the grid data
        print("Found z_average.vtmb, loading it...\n")
        dataset = pv.read( filename )
        dataset = dataset.cell_data_to_point_data().combine()
        print(dataset.array_names)
        print("Finished data preparation!")
    
    return dataset, x_pfmax
    
def pv_visualize( dataset, varname, clipbox, cbar_ticks, x_pfmax ):
    
    pl = pv.Plotter(off_screen=True, window_size=[1920,1080], border=False)
    
    cmap    = plt.get_cmap('coolwarm',51)
    clim    = [cbar_ticks[0], cbar_ticks[-1]]
    
    sepline = dataset.contour( [0.0], scalars='u' )
    sonline = dataset.contour( [1.0], scalars='mach' )
    pl.add_mesh(dataset, scalars=varname, show_scalar_bar=False, cmap=cmap,
                clim=clim)
    if sepline.n_points > 0:
        pl.add_mesh(sepline, color='yellow', line_width=2.0 )
    pl.add_mesh(sonline, color='lime',  line_width=2.0 )
    pl.view_vector([0.0,0.0,1.0],viewup=[0.0,1.0,0.0])
    pl.camera.tight()
    image = crop_to_rect_map(pl.screenshot(return_img=True), buff=100)
    pl.close()
    
    fig, ax = plt.subplots(figsize=(12.8,7.2))
    img = ax.imshow(image, extent=clipbox[:4], cmap=cmap, clim=clim)
    
    ax.plot( x_pfmax, 0.0, '*', color='cyan' , markersize=20 )

    lines = list()
    x_sep = min( sepline.points[:,0] )
    x_att = max( sepline.points[:,0] )
    bbox = [x_sep, clipbox[1], clipbox[2], clipbox[3]]
    line = dividing_streamline( dataset, np.array([-1.0,0.3,0.0]), bbox,step=0.002, forward=True )
    lines.append( line )
    line = dividing_streamline( dataset, np.array([x_sep,0.005,0.0]), bbox, step=0.002, forward=True )
    lines.append( line )
    
    with open('dividing_streamline.pkl', 'wb') as f:
        pickle.dump([line], f)   # only the last line is saved!
    
    if lines is not None:
        for line in lines:
            ax.plot( line[:,0], line[:,1], color='black', linewidth=1.0 )
    
    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')

    ax.set_xlim(clipbox[0], clipbox[1])
    ax.set_ylim(clipbox[2], 6.0)
    ax.set_aspect('equal')

    cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5, extend='both' ) 
    cbar.ax.set_ylabel( varname, loc='center', labelpad=30)
    cbar.ax.set_xticks(cbar_ticks)
    
    plt.show()
    

def test_fetch_uv():
    """
    test the fetch_uv function
    """
    
    case_dir = '/home/wencan/temp/220927/' 
    dirs     = Directories( case_dir )
    
    dataset  = pv.read( dirs.pp_z_average + '/z_average.vtmb' )
    dataset  = dataset.cell_data_to_point_data().combine()
    
    point_start = np.array([-3.57135,1.007,0.0])
    
    u,v = fetch_uv(dataset, point_start)
    
    print(f"u: {u}, v: {v}")


# =============================================================================
if __name__ == "__main__":

    main()

