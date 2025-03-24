#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   bubble_size_from_pdf.py
@Time    :   2024/05/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Script of computing the bubble size from the separation pdf.

             Need inca_grid.bin, cutcells_setup.dat, snapshot_container.pkl and wall_dist/ in the same directory.
'''


import os
import sys
import pickle
import numpy             as     np
import pandas            as     pd
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid         import GridData
from   vista.timer        import timer
from   vista.params       import Params
from   vista.snapshot     import Snapshot
from   vista.statistic    import StatisticData
from   vista.directories  import Directories
from   vista.tools        import get_filelist
from   vista.tools        import lin_grow
from   vista.plane_analy  import pv_interpolate
from   vista.plot_setting import set_plt_rcparams

set_plt_rcparams(latex=False,fontsize=15)

def main():

# =============================================================================

    case_dir = '/home/wencan/temp/smooth_adiabatic/'
    var_out  = ['u','n_sep','pdf_sep','wd']
    rescale  = [-50.4, 0.0, 0.0, 5.2, 5.2, 5.2]

# =============================================================================

    dirs      = Directories( case_dir )
    params    = Params( dirs.case_para_file )
    roughwall = params.roughwall

    # load the grid data

    grd       = GridData( dirs.grid )
    grd.read_grid()
    grd.cell_volume()

    # load the snapshot container of separation PDF

    os.chdir( dirs.pp_bubble )
    with open(f"{dirs.pp_bubble}/snapshot_container.pkl",'rb') as f:
        snapshot_container: Snapshot = pickle.load(f)

    blocklist = snapshot_container.bl_nums

    # load the statistic data

    stat      = StatisticData( dirs.statistics )
    stat.read_statistic( block_list=blocklist, vars_in=['u'] )

    # load the wall distance field if roughwall 

    if roughwall:
        wd_snap_file = get_filelist( dirs.wall_dist, key='snapshot.bin')[0]
        with timer("load wall distance snapshot data"):
            wd_snap = Snapshot( wd_snap_file )
            wd_snap.read_snapshot( var_read=['wd'] )
        
        snapshot_container.copy_var_from( wd_snap, ['wd'] )

    else:
        
        for bl in snapshot_container.snap_data:
            bl_num       = bl.num
            g            = grd.g[bl_num-1]
            X, Y, Z      = np.meshgrid( g.gx, g.gy, g.gz, indexing='ij' )
            bl.df['wd']  = Y.T.flatten()

    # copy statistic u to snapshot container

    for bl in snapshot_container.snap_data:
        bl_num       = bl.num
        statbl       = stat.bl[stat.bl_nums.index(bl_num)]
        bl.df['u']   = statbl.df['u']

    # set cutcells data volume fraction

    if roughwall:
        cc_df = pd.read_csv( dirs.cc_setup, delimiter=r'\s+' )
        cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], 
                    inplace=True)
    else: cc_df = None


    # compute the bubble size by the definition of 50% pdf

    with timer("compute bubble size"):
        
        sep1,sep2 = snapshot_container.compute_bubble_volume_pdf( grd, cc_df, 
                                   roughwall=roughwall, opt=1, y_threshold=0.0 )
        sep3,sep4 = snapshot_container.compute_bubble_volume_pdf( grd, cc_df, 
                                   roughwall=roughwall, opt=2, y_threshold=0.0 )

        df_bubble = pd.read_csv('bubble_size.dat', delimiter=r'\s+')
        sep5      = np.mean( df_bubble['bubble_volume'] )
        sep6      = np.mean( df_bubble['bubble_volume_thr'] )

        print(f"bubble size by 50% pdf   (threshold y>0): {sep1:.2f} ({sep2:.2f})")
        print(f"bubble size PDF          (threshold y>0): {sep3:.2f} ({sep4:.2f})")
        print(f"bubble size time average (threshold y>0): {sep5:.2f} ({sep6:.2f})")
        
        with open('bubble_size_pdf.dat','w') as f:
            f.write(f"bubble size by 50% pdf  (threshold y>0): {sep1:.2f} ({sep2:.2f})\n")
            f.write(f"bubble size PDF         (threshold y>0): {sep3:.2f} ({sep4:.2f})\n")
            f.write(f"bubble size time average(threshold y>0): {sep5:.2f} ({sep6:.2f})\n")

# =============================================================================
# --- show bubble outline computed from p.d.f. and mean u

    snapshot_container.grid3d = grd
    dataset = snapshot_container.spanwise_average( blocklist, ['u','pdf_sep'], rescale=rescale )
    dataset = pv.MultiBlock( dataset )

    post_process_dataset( dataset, ['u','pdf_sep'] )

    # dataset = dataset.cell_data_to_point_data().combine()

    # p = pv.Plotter()
    # sep1 = dataset.contour( [0.5],scalars='pdf_sep')
    # sep2 = dataset.contour( [0.0],scalars='u')
    # p.add_mesh( sep1, color='yellow', line_width=4.0 )
    # p.add_mesh( sep2, color='green',  line_width=3.0 )
    # p.show()
        
    # write the PDF file into tecplot szplt format.

    # snapshot_container.grid3d = grd
    # snapshot_container.write_vtm( 'pdf_sep.vtm', vars=var_out, buff=2 )

    #snapshot_container.write_szplt( 'pdf_sep.szplt', buff=2 )


# =============================================================================
# --- post process the dataset to show the bubble outline

def post_process_dataset( dataset:pv.MultiBlock, vars_out ):

    px   = np.linspace( -15, 10, 301, endpoint=True )
    py,_ = lin_grow( 0.0, 0.03, 1.04, upbound=6.0 )
    
    # points near interface where mesh density jumps will lead to holes
    
    # px   = np.linspace( -4.2, -3.4, 11,  endpoint=True ) 
    # py   = np.linspace(  3.3,  3.5, 101, endpoint=True )
    pz   = np.array([0.0])

    dataset = dataset.cell_data_to_point_data().combine()

    df = pv_interpolate( dataset, vars_out,[px,py,pz])

    x        = np.array( df['x']      ).reshape( (len(py),len(px)) )
    y        = np.array( df['y']      ).reshape( (len(py),len(px)) )
    u        = np.array( df['u']      ).reshape( (len(py),len(px)) )
    pdf_sep  = np.array( df['pdf_sep']).reshape( (len(py),len(px)) )

    fig, ax = plt.subplots(1,1,figsize=(12.8,7.2))
    c       = ax.contourf( x, y, pdf_sep,  levels=np.linspace(0,1.0,51),  cmap='Reds', extend='both')
    csepu   = ax.contour(  x, y, u,        linestyles='solid', levels=[0.0],  
                           colors='blue',  linewidths=1.5,     zorder=10)
    cseppdf = ax.contour(  x, y, pdf_sep,  linestyles='dashed',levels=[0.5],  
                           colors='lime',  linewidths=1.5,     zorder=10)
    
    ax.set_aspect('equal')
    
    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')
    
    ax.spines[:].set_linewidth(1.0)
    ax.tick_params( direction='out', length=5, width=1.5)
    ax.set_xlim([-12,5])
    ax.set_ylim([0,3])
    
    cbar = fig.colorbar( c, 
                         ax=ax, 
                         pad=0.20,
                         shrink=0.5,
                         orientation='horizontal',
                         ticks=np.linspace(0,1.0,5))
    cbar.ax.tick_params( direction='in',
                         length=5,
                         width=1.5)
    cbar.ax.set_ylabel('p.d.f.', loc='center', labelpad=20)
    cbar.outline.set_linewidth(1.5)
    
    plt.savefig('pdf_sep_averaged.png', dpi=300, bbox_inches='tight')
    plt.show()

    plt.close()
    

# =============================================================================

if __name__ == '__main__':
    main()
