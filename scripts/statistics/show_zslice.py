#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_zslice.py
@Time    :   2025/03/28 
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
    case_dir   = '/home/wencan/temp/241030/'
    vars_read  = ['u','v','w','p','T','pp','uu','vv','ww','uv','uw','vw']
    vars_out   = ['u','v','w','p','T','mach','grad_rho','DS','tke','uu',
                 'u`u`','v`v`','w`w`','u`v`','p`']
    rescale    = [-50.4, 0.0, 0.0, 5.2, 5.2, 5.2]
    clipbox    = [-11.5, -7.5, 0, 0.5, -1, 1]
    cbar_ticks = np.linspace(0.0,200000,4, endpoint=True)
    
    dataset, x_pfmax = data_preparation( case_dir, vars_read, vars_out, rescale )
    dataset = dataset.clip_box( clipbox, invert=False )
    
    pv_visualize( dataset, 'uu', clipbox, cbar_ticks, x_pfmax )

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
    stat2d.compute_vars(      blocklist, ['mach','RS','p`'] )
    stat2d.compute_gradients( blocklist, ['grad_rho','DS'] )
    dataset = pv.MultiBlock(stat2d.create_vtk_multiblock(blocklist, vars_out,rescale=rescale))
    dataset = dataset.cell_data_to_point_data().combine()
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
        pl.add_mesh(sepline, color='yellow', line_width=4.0 )
    pl.add_mesh(sonline, color='lime',  line_width=4.0 )
    pl.view_vector([0.0,0.0,1.0],viewup=[0.0,1.0,0.0])
    pl.camera.tight()
    image = crop_to_rect_map(pl.screenshot(return_img=True), buff=100)
    pl.close()
    
    fig, ax = plt.subplots(figsize=(12.8,7.2))
    img = ax.imshow(image, extent=clipbox[:4], cmap=cmap, clim=clim)
    
    ax.plot( x_pfmax, 0.0, '*', color='cyan' , markersize=20 )
    
    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')

    cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5, extend='both' ) 
    cbar.ax.set_ylabel( varname, loc='center', labelpad=30)
    cbar.ax.set_xticks(cbar_ticks)
    
    plt.show()
    


# =============================================================================
if __name__ == "__main__":

    main()

