#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   boundary_edge_streamline.py
@Time    :   2025/05/07 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   show and save the boundary edge streamline
'''

import os
import sys
import pickle
import numpy             as     np
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.params      import Params
from   vista.directories import Directories
from   vista.tools       import crop_to_rect_map
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams( fontsize=20 )

def main():
    case_dir   = '/home/wencan/temp/smooth_adiabatic/' 
    clipbox    = [-15, 10, 0, 10.0, -1, 1]
    cbar_ticks = np.linspace(0.0,2.0,5, endpoint=True)
    
    dataset, x_pfmax = data_preparation( case_dir )
    dataset = dataset.clip_box( clipbox, invert=False )
    
    dirs    = Directories( case_dir )
    os.chdir( dirs.pp_z_average )
    
    pv_visualize( dataset, 'mach', clipbox, cbar_ticks, x_pfmax )

    
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
    for h in np.linspace(2.0,2.0,1, endpoint=True):
        line = boundary_edge_streamline( dataset, np.array([-14.5,h,0.0]), step=0.05, forward=True )
        lines.append( line )
    
    with open('boundary_edge_streamline.pkl', 'wb') as f:
        pickle.dump([line], f)
    
    if lines is not None:
        for line in lines:
            ax.plot( line[:,0], line[:,1], color='black', linewidth=1.0 )
    
    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')

    ax.set_xlim(clipbox[0], clipbox[1])
    ax.set_ylim(clipbox[2], clipbox[3])
    ax.set_aspect('equal')

    cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5, extend='both' ) 
    cbar.ax.set_ylabel( varname, loc='center', labelpad=30)
    cbar.ax.set_xticks(cbar_ticks)
    
    plt.show()
    


def boundary_edge_streamline( dataset:pv.UnstructuredGrid, point_start:np.ndarray, step=0.1, forward=True):

    """
    compute the dividing streamline
    """

    u,v = fetch_uv(dataset, point_start)
    pts = [point_start]
    
    i = 0
    while True:
        
        point_start = step_streamline( point_start, u, v, step=step, forward=forward )
        pts.append( point_start )
        u, v = fetch_uv(dataset, point_start)
        i += 1
        
        if i > 3000:
            break
        
        if point_start[0] < -15.0 or point_start[0] > -4.0 or \
           point_start[1] < 0.0   or point_start[1] > 10.0:
            break
        
    return np.array(pts)

def step_streamline( point_start, u, v, step=0.05, forward=True):
    """
    step the streamline
    """

    if not forward:
        step = -step
    
    v_mag = np.sqrt(u**2 + v**2)
    
    x = point_start[0] + step*u/v_mag
    y = point_start[1] + step*v/v_mag
    z = point_start[2]
    
    end_point = np.array([x,y,z])
    
    return end_point
    

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


# =============================================================================
if __name__ == "__main__":

    main()

