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
import numpy             as     np
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.params      import Params
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.tools       import crop_to_rect_map
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams( fontsize=20 )

def main():
    case_dir   = '/home/wencan/temp/smooth_adiabatic/'
    vars_read  = ['u','v','w','p','T']
    vars_out   = ['u','v','w','p','T','mach','a','DS']
    rescale    = [-50.4, 0.0, 0.0, 5.2, 5.2, 5.2]
    clipbox    = [-15, 10, 0, 10.0, -1, 1]
    cbar_ticks = np.linspace(0.0,0.8,5, endpoint=True)
    
    dataset, x_pfmax = data_preparation( case_dir, vars_read, vars_out, rescale )
    dataset = dataset.clip_box( clipbox, invert=False )
    
    lines = list()
    
    for x in np.linspace( -7.0, -6.0, 12, endpoint=True ):
    
        line1 = c_line(dataset, np.array([x,2.0, 0.0]), step=0.11, forward=False, left=True)
        line2 = c_line(dataset, np.array([x,2.0, 0.0]), step=0.11, forward=False, left=False)
        lines.append( line1 )
        lines.append( line2 )
    
    pv_visualize( dataset, 'DS', clipbox, cbar_ticks, x_pfmax, lines )

    
def c_line(dataset:pv.UnstructuredGrid, point_start:np.ndarray, step=0.1, forward=False, left=True):
    """
    compute the reversed characteristic line
    """
    u, v, a = fetch_uva(dataset, point_start)
    pts = [point_start]
    
    i = 0
    while True:
        
        point_start = step_c(point_start, u, v, a, step=step, forward=forward, left=left)
        pts.append(point_start)
        u, v, a     = fetch_uva(dataset, point_start)
        i += 1

        if a > np.sqrt(u**2+v**2):
            print("reach the subsonic region!")
            break
        
        if i > 100:
            break
        
        if point_start[0] < -15 or point_start[0] > 10 or \
           point_start[1] < 0.0 or point_start[1] > 10.0:
            break
        
        print(i, point_start, u, v, a)
    
    return pts
    
    
    
def step_c( point_start, u, v, a, step=0.1, forward=True, left=True ):
    """
    step forward/backward in space
    """
    # - compute the backward step
    theta = np.arctan(v/u)
    mu    = np.arcsin(a/np.sqrt(u**2+v**2))
    
    if left:
        beta = theta+mu
    else:
        beta = theta-mu
    
    if not forward:
        step = -step
    
    dx = step * np.cos(beta)
    dy = step * np.sin(beta)
    
    point_end = np.array([point_start[0]+dx, point_start[1]+dy, point_start[2]])

    return point_end    
    

def fetch_uva( dataset:pv.UnstructuredGrid, point:np.ndarray ):
    """
    fetch the u, v, a at the point
    """
    pt=pv.PolyData(point)
    results = pt.sample( dataset )
    
    return results['u'][0], results['v'][0], results['a'][0]
    
    
def pv_visualize( dataset, varname, clipbox, cbar_ticks, x_pfmax, lines ):
    
    pl = pv.Plotter(off_screen=True, window_size=[1920,1080], border=False)
    
    cmap    = plt.get_cmap('coolwarm',51)
    clim    = [cbar_ticks[0], cbar_ticks[-1]]
    
    sepline = dataset.contour( [0.0], scalars='u' )
    sonline = dataset.contour( [1.0], scalars='mach' )
    pl.add_mesh(dataset, scalars=varname, show_scalar_bar=False, cmap=cmap,
                clim=clim)
    if sepline.n_points > 0:
        pl.add_mesh(sepline, color='yellow', line_width=4.0 )
    pl.add_mesh(sonline, color='lime',  line_width=4.0 )
    pl.view_vector([0.0,0.0,1.0],viewup=[0.0,1.0,0.0])
    pl.camera.tight()
    image = crop_to_rect_map(pl.screenshot(return_img=True), buff=100)
    pl.close()
    
    fig, ax = plt.subplots(figsize=(12.8,7.2))
    img = ax.imshow(image, extent=clipbox[:4], cmap=cmap, clim=clim)
    
    ax.plot( x_pfmax, 0.0, '*', color='cyan' , markersize=20 )
    
    for line in lines:
        line = np.array(line)
        ax.plot( line[:,0], line[:,1], color='teal', linewidth=1.0 )
    
    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')
    
    ax.set_xlim(clipbox[0], clipbox[1])
    ax.set_ylim(clipbox[2], clipbox[3])
    ax.set_aspect('equal')

    cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5, extend='both' ) 
    cbar.ax.set_ylabel( varname, loc='center', labelpad=30)
    cbar.ax.set_xticks(cbar_ticks)
    
    plt.show()
    

def data_preparation( case_dir, vars_read, vars_out, rescale ):

    dirs     = Directories( case_dir )
    
    params   = Params( dirs.case_para_file )
    x_pfmax  = params.x_pfmax

    # - check if dirs.stat_zslice exists

    if not os.path.exists( dirs.stat_zslice ):
        print(f"Directory {dirs.stat_zslice} does not exist.")
        print("Please run 'scripts/statistics/write_zslice.py' first!")
        sys.exit()
    
    # - read in the grid data
    grid = GridData( dirs.grid )
    grid.read_grid()
    blocklist,_= grid.select_sliced_blockgrids( 'Z', 0.001 )
    
    # - read in the z-slice statistics

    stat2d = StatisticData( dirs.stat_zslice )
    stat2d.grid3d = grid
    stat2d.read_statistic(    blocklist, vars_in=vars_read )
    stat2d.compute_vars( blocklist, ['mach'] )
    stat2d.compute_gradients( blocklist, ['grad_rho','DS'] )
    dataset = pv.MultiBlock(stat2d.create_vtk_multiblock(blocklist, vars_out,buff=2,rescale=rescale))
    dataset = dataset.cell_data_to_point_data().combine()
    print("Finished data preparation!")
    
    return dataset, x_pfmax

# =============================================================================
if __name__ == "__main__":

    main()

