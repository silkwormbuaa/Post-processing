#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_zslice.py
@Time    :   2025/04/10
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Show variable in the incipient interaction plane
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
    case_dirs  = ['/home/wencan/temp/smooth_mid/']
    vars_read  = ['u','v','w','p','T','pp','uu','vv','ww','uv','uw','vw']
    vars_out   = ['u','v','w','p','T','mach','grad_rho','DS','tke',
                 'u`u`','v`v`','w`w`','u`v`','p`']
    clipbox    = [-100, 100, 0.0, 1.0, 0, 0.5]

    for case_dir in case_dirs:
        print(f"Start processing {case_dir}.\n")
        
        dirs     = Directories( case_dir )
        
        dataset  = pv.read( dirs.stat_incip + '/stat_incip.vtm' )
        dataset  = dataset.cell_data_to_point_data().combine()
        dataset  = dataset.clip_box( clipbox, invert=False )
        
        cmap = 'coolwarm'
        cbar_ticks = np.linspace(0,1500, 5, endpoint=True)
        pv_visualize( dataset, 'p`', clipbox, cbar_ticks, cmap )
        
        cmap = 'coolwarm'
        cbar_ticks = np.linspace(0,8000, 4, endpoint=True)
        pv_visualize( dataset, 'tke', clipbox, cbar_ticks, cmap )

        cmap = 'coolwarm_r'
        cbar_ticks = np.linspace(-800,0, 4, endpoint=True)
        pv_visualize( dataset, 'u`v`', clipbox, cbar_ticks, cmap )

    
def pv_visualize( dataset, varname, clipbox, cbar_ticks, c_map ):
    
    pl = pv.Plotter(off_screen=True, window_size=[1920,1080], border=False)
    
    cmap    = plt.get_cmap(c_map,51)
    clim    = [cbar_ticks[0], cbar_ticks[-1]]
    
    sepline = dataset.contour( [0.0], scalars='u' )
    sonline = dataset.contour( [1.0], scalars='mach' )
    pl.add_mesh(dataset, scalars=varname, show_scalar_bar=False, cmap=cmap,
                clim=clim)
    if sepline.n_points > 0:
        pl.add_mesh(sepline, color='yellow', line_width=4.0 )
    pl.add_mesh(sonline, color='lime',  line_width=4.0 )
    pl.view_vector([1.0,0.0,0.0],viewup=[0.0,1.0,0.0])
#    pl.camera.tight()
#    pl.add_axes()
    image = crop_to_rect_map(pl.screenshot(return_img=True))
    pl.close()
    
    fig, ax = plt.subplots(figsize=(12.8,7.2))
    img = ax.imshow(image, extent=clipbox[4:]+clipbox[2:4], cmap=cmap, clim=clim)
    
    ax.set_xlabel(r'$z/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')
    ax.set_aspect('equal')

    cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5, extend='both' ) 
    cbar.ax.set_ylabel( varname, loc='center', labelpad=30)
    cbar.ax.set_xticks(cbar_ticks)
    
    plt.show()
    plt.close()
    


# =============================================================================
if __name__ == "__main__":

    main()

