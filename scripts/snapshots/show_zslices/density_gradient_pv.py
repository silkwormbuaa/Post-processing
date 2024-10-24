#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   density_gradient_pv.py
@Time    :   2024/08/21 
@Author  :   Wencan WU 
@Version :   2.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   read in a snapshot and show it with pyvista and matplotlib
'''

off_screen = False

if off_screen:
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb(width=1920, height=1080)
    vdisplay.start()
    
import os
import sys
import time
import pyvista            as     pv
import matplotlib.pyplot  as     plt
from   mpi4py             import MPI

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid         import GridData
from   vista.timer        import timer
from   vista.snapshot     import Snapshot
from   vista.directories  import Directories
from   vista.tools        import crop_border
from   vista.tools        import get_filelist
from   vista.tools        import distribute_mpi_work
from   vista.directories  import create_folder
from   vista.plot_setting import cpos_callback
from   vista.plot_setting import set_plt_rcparams
#from   vista.log         import Logger
#sys.stdout = Logger( os.path.basename(__file__) )

# - build MPI communication environment

comm    = MPI.COMM_WORLD
rank    = comm.Get_rank()
n_procs = comm.Get_size()

# =============================================================================

case_dir = '/home/wencan/temp/231124'
out_dir  = '/home/wencan/temp/231124/postprocess/snapshotZ/'
clipbox  = [-20,12,0,10,-1,1]

# =============================================================================

# =============================================================================
# preparation
# =============================================================================

dirs      = Directories( case_dir )
grid3d    = GridData()
snapfiles = None

if rank == 0:

    create_folder( out_dir )
    
    snapfiles = get_filelist( dirs.snp_dir, 'snapshot_Z_001.bin' )
    print(f"I am root, just found {len(snapfiles)} snapshot Z 001 files.")
    
    grid3d = GridData( dirs.grid )
    grid3d.read_grid()

snapfiles  = comm.bcast( snapfiles, root=0 )
grid3d     = comm.bcast( grid3d,    root=0 )

n_snaps    = len(snapfiles)
i_s, i_e   = distribute_mpi_work( n_snaps, n_procs, rank )
snapfiles  = snapfiles[i_s:i_e]

print(f"I am processor {rank:05d}, I take care of {len(snapfiles):5d} snapshots.")
sys.stdout.flush()

vars       = ['u','v','w','p','T','grad_rho','grad_p']
rescale    = [-50.4,.0,.0,5.2,5.2,5.2]

os.chdir( out_dir )
clock = timer("show snapshot z slices")

# -- loop over the snapshots

for i, snap_file in enumerate(snapfiles):

    snap        = Snapshot( snap_file )
    snap.grid3d = grid3d
    snap.read_snapshot()
    snap.compute_gradients( grads=['grad_rho','grad_p'] )
    
# ---- pass snapshot to pyvista

    dataset = pv.MultiBlock( snap.create_vtk_multiblock( vars, mode='symmetry', rescale=rescale ) )

    p = pv.Plotter( off_screen=True, window_size=[1920,1080], border=False )

    dataset = dataset.cell_data_to_point_data().combine()
    dataset = dataset.clip_box(clipbox, invert=False)
    dataset.set_active_scalars('grad_rho')

    cmap    = plt.get_cmap('viridis',51)
    cmap.set_over('red')
    cmap.set_under('blue')

    p.add_mesh( dataset, 
                cmap=cmap, 
                clim=[0.0,0.8],
                show_scalar_bar=False )

    p.view_vector([0.0,0.0,1.0],viewup=[0.0,1.0,0.0])
    p.camera.tight()
    image = p.screenshot(return_img=True)
    p.close()

    # >>> show bounds with pyvista: ugly font and hard to control the font size,
    #                               but useful in interactive mode
    # =============================================================================
    #
    # p.update_bounds_axes()
    # p.show_bounds( dataset, 
    #                bounds=[-20,12,0,10,0,0],
    #                ticks='both', show_zaxis=False,
    #                xtitle=r'$(x-x_{imp})/\delta_0$',
    #                ytitle=r'$y/\delta_0$',
    #                location='front',
    #                font_size=15)
    # p.camera_position = [(-4,5,40),(-4,5,0),(0,1,0)]
    # p.add_title(f"t = {snap.itime:.2f} s", font_size=15, color='white')
    # p.show_axes()
    # p.show()
    #
    # =============================================================================

    image = crop_border(image)
    set_plt_rcparams(fontsize=15)

    fig, ax = plt.subplots(figsize=(12.8,7.2))

    img = ax.imshow(image, extent=clipbox[:4], cmap=cmap, clim=[0.0,0.8])

    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')

    cbar = fig.colorbar( img, orientation='horizontal', ax=ax, shrink=0.5 ) 
    cbar.ax.set_ylabel( r'$\nabla \rho$', loc='center')

    plt.title(f"t = {snap.itime:.2f} ms", loc='center',y=1.05)

    if off_screen:
        plt.savefig(f"grad_rho_{i_s+i:06d}.png", dpi=150)
    else:
        plt.show()
    
    plt.close()

# -- print the progress

    progress = (i+1)/len(snapfiles)
    print(f"Rank:{rank:05d},{i+1}/{len(snapfiles)} is done. " + clock.remainder(progress))
    print("------------------\n")
    sys.stdout.flush()

comm.barrier()

if rank == 0:
    
    print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    sys.stdout.flush() 
