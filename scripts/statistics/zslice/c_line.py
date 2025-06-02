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

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.params      import Params
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.tools       import crop_to_rect_map
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams( fontsize=30 )

def main():
    case_dir   = '/home/wencan/temp/smooth_mid/'
    vars_read  = ['u','v','w','p','T']
    vars_out   = ['u','v','w','p','T','mach','a','DS']
    rescale    = [-50.4, 0.0, 0.0, 5.2, 5.2, 5.2]
    clipbox    = [-12.0, 5.0, 0, 6.0, -1, 1]
    cbar_ticks = np.linspace(0,0.8,5, endpoint=True)
    
    dataset, x_pfmax = data_preparation( case_dir, vars_read, vars_out, rescale )
    dataset = dataset.clip_box( clipbox, invert=False )
    
    lines = list()
    
    # smooth_mid : -6.9, -6.4
    # 231124     : -9.9, -9.4
    # 241030     : -8.4, -7.9
    
    for x in np.linspace( -6.9, -6.5, 12, endpoint=True ):
    
        line1 = c_line(dataset, np.array([x,2.0, 0.0]), step=0.05, forward=False, left=True)
        line2 = c_line(dataset, np.array([x,2.0, 0.0]), step=0.11, forward=False, left=False)
        lines.append( line1 )
#        lines.append( line2 )
    
    dirs   = Directories( case_dir )
    params = Params( dirs.case_para_file )
    
    x_start = -2.5
    x_stop  = params.x_pfmax
    line_wall = c_line_wall(dataset, np.array([x_start,0.05, 0.0]), x_stop=x_stop, step=0.05, forward=False)
    
    pv_visualize( dataset, 'DS', clipbox, cbar_ticks, x_pfmax, lines, case_dir )

    
def c_line(dataset:pv.UnstructuredGrid, point_start:np.ndarray, step=0.1, forward=False, left=True):
    """
    compute the reversed characteristic line
    """
    u, v, a = fetch_uva(dataset, point_start)
    pts = [point_start]
    time= 0.0
    
    i = 0
    while True:
        
        point_start = step_c(point_start, u, v, a, step=step, forward=forward, left=left)
        pts.append(point_start)
        
        ds = np.sqrt( (point_start[0]-pts[-2][0])**2 + (point_start[1]-pts[-2][1])**2 )
        time += ds / np.sqrt(u**2+v**2+a**2)
        
        u, v, a     = fetch_uva(dataset, point_start)
        i += 1

        if a > np.sqrt(u**2+v**2):
            print("reach the subsonic region!")
            break
        
        if i > 1000:
            break
        
        if point_start[0] < -15 or point_start[0] > 10 or \
           point_start[1] < 0.0 or point_start[1] > 10.0:
            break
        
#        print(i, point_start, u, v, a)
    
    print(f"total time: {time*5.2:.4f} ms")
    
    return pts
    
def c_line_wall(dataset:pv.UnstructuredGrid, point_start:np.ndarray, x_stop:float, step=0.1, forward=False):
    """
    compute the reversed characteristic line
    """
    u, v, a = fetch_uva(dataset, point_start)
    pts = [point_start]
    time= 0.0
    
    i = 0
    while True:
        
        point_start = step_c_wall(point_start, step=step, forward=forward)
        pts.append(point_start)
        
        ds = step
        if forward:
            time += ds / (a+u)
        else:
            time += ds / (a-u)
            print(f"{ds / (a-u):.7f} ms")
        
        u, v, a     = fetch_uva(dataset, point_start)
        i += 1
        
        if i > 1000:
            break
        
        if point_start[0] < x_stop:
            break
        
#        print(i, point_start, u, v, a)
    
    print(f"total time on wall: {time*5.2:.4f} ms")
    
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


def step_c_wall( point_start, step=0.1, forward=True ):
    """
    step forward/backward in space
    """
    # - compute the backward step
    
    
    if not forward:
        step = -step
    
    dx = step 
    
    point_end = np.array([point_start[0]+dx, point_start[1], point_start[2]])

    return point_end    


def fetch_uva( dataset:pv.UnstructuredGrid, point:np.ndarray ):
    """
    fetch the u, v, a at the point
    """
    pt=pv.PolyData(point)
    results = pt.sample( dataset )
    
    return results['u'][0], results['v'][0], results['a'][0]
    
    
def pv_visualize( dataset, varname, clipbox, cbar_ticks, x_pfmax, clines, dirs ):
    
    dirs = Directories( dirs )
    
    pl = pv.Plotter(off_screen=True, window_size=[1920,1080], border=False)
    
    cmap    = plt.get_cmap('Greys_r',51)
    clim    = [cbar_ticks[0], cbar_ticks[-1]]
    
    sepline = dataset.contour( [-0.1], scalars='u' )
    sonline = dataset.contour( [1.0] , scalars='mach' )
    pl.add_mesh(dataset, scalars=varname, show_scalar_bar=False, cmap=cmap,
                clim=clim,lighting=False)
    # pl.add_mesh(sonline, color='lime',  line_width=4.0 )
    # if sepline.n_points > 0:
    #     pl.add_mesh(sepline, color='red', line_width=4.0 )
    pl.view_vector([0.0,0.0,1.0],viewup=[0.0,1.0,0.0])
    pl.camera.tight()
    image = crop_to_rect_map(pl.screenshot(return_img=True), buff=100)
    pl.close()
    
    fig, ax = plt.subplots(figsize=(12.8,7.2))
    img = ax.imshow(image, extent=clipbox[:4], cmap=cmap, clim=clim)
    
    ax.plot( x_pfmax, 0.0, '*', color='black' , markersize=20, zorder=10 )
    
    for i, line in enumerate(clines):
        line = np.array(line)
        if i == 5:
            ax.plot( line[:,0], line[:,1], color='orange', ls='--',linewidth=2.0, zorder=9 )
        else: ax.plot( line[:,0], line[:,1], color='teal', ls='--',linewidth=1.5 )
    
    sonlines = pickle.load( open( dirs.pp_z_average + '/sonicline.pkl', 'rb' ) )
    for line in sonlines:
        line = np.array(line)
        ax.plot( line[:,0], line[:,1], color='lime', linewidth=2.0 )
    
    seplines = pickle.load( open( dirs.pp_z_average + '/sepline.pkl', 'rb' ) )
    for line in seplines:
        line = np.array(line)
        ax.plot( line[:,0], line[:,1], color='red', linewidth=2.0 )
    
    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')
    
    ax.set_xlim(clipbox[0], clipbox[1])
    ax.set_ylim(clipbox[2], clipbox[3])
    ax.set_aspect('equal')
    ax.spines[:].set_linewidth(2.0)
    ax.tick_params( direction='out', length=10, width=2.0)

    cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5, extend='both', pad=0.2 ) 
    cbar.ax.set_ylabel( varname, loc='center', labelpad=50)
    cbar.ax.set_xticks(cbar_ticks)
    # move cbar away from the axes of image
    cbar.ax.tick_params( direction='in', length=5, width=1.5 )
    cbar.outline.set_linewidth(1.5)
    plt.tight_layout()
    
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

