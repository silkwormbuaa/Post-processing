#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_bubble_pdf.py
@Time    :   2025/03/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

off_screen = False

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()

import os
import sys
import pyvista            as     pv
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.snapshot    import Snapshot
from   vista.directories import Directories
from   vista.tools       import crop_border
from   vista.tools       import get_filelist
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(latex=False,fontsize=15)


def main():
   
    case_dir = '/home/wencan/temp/250120/'
    clipbox  = [-32.8,102.4,-2,15,-11,11]
    dirs = Directories( case_dir )
    
    show_bubble_pdf( dirs, clipbox )
    

def show_bubble_pdf( dirs:Directories, clipbox:list ):

    dataset   = pv.read( dirs.pp_bubble + '/pdf_sep.vtm' )
    print( dataset )
    pointdata = dataset.cell_data_to_point_data().combine()
    print( pointdata )
    dataset = dataset.clip_box( clipbox , invert=True )

    p = pv.Plotter(off_screen=True, window_size=[1920,1080], border=False)
    
    wall   = pointdata.contour( [0.01], scalars='wd' )
    bubble = pointdata.contour( [0.5], scalars='pdf_sep' )
    
    p.add_mesh( wall, color='grey' )
    p.add_mesh( bubble, color='blue' )
    
    p.view_vector( [0.0, 1.0, 0.0], viewup=[0.0,0.0,-1.0])
    p.enable_parallel_projection()

    # p.camera.tight(view='xz')
    # p.show()
    
    print( "screenshot..." )
    
    image = p.screenshot( return_img=True )
    p.close()
    
    image = crop_border( image )
    
    fig, ax = plt.subplots(figsize=(12.8,7.2))
    
    img = ax.imshow( image, extent=[-16,10,-2,2] )
    
    
    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')
    
    plt.show()
    
    plt.close()
    
    


# =============================================================================
if __name__ == "__main__":

    main()

