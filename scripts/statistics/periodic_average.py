#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   periodic_average.py
@Time    :   2025/02/25 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import numpy             as     np
import pandas            as     pd
import pyvista           as     pv
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.params      import Params
from   vista.tools       import get_filelist
from   vista.tools       import define_wall_shape
from   vista.directories import create_folder
from   vista.plane_analy import pv_interpolate

plt.rcParams["text.usetex"]         = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
plt.rcParams['font.family']         = "Times New Roman"
plt.rcParams['font.size']           = 30

def main():
# =============================================================================

    case_folder = '/home/wencan/temp/smooth_adiabatic/'
    
    stat_file   = 'stat_xslice_-12.bin'
    outfolder   = 'yz_planes_-12'
    loc         = -12.0
    
    loc         = loc * 5.2 + 50.4
    
    dirs        = Directories( case_folder )
    params      = Params( dirs.case_para_file )
    
    vars_read   = ['u','v','w','p','T','uu','vv','ww','uv','pp']
    vars_output = ['u',    'v',   'w',    'T',
                   'mach', 'tke', 'u`u`', 'u`v`',
                   'p`']
    cbar_range  = [[0.0,1.0],[-4.0,4.0],[-2.0,2.0], [1.0,1.8],
                   [0,2.0],   [0,4.0],   [0,3.0],   [-3.6,0],
                   [0,0.08]]
    cbar_label  = [r'$\langle u \rangle / u_{\infty}$',
                   r'$\langle v \rangle / u_{\infty}\cdot 100$',
                   r'$\langle w \rangle / u_{\infty}\cdot 100$',
                   r'$\langle T \rangle / T_{\infty}$',
                   r'$\mathrm{Mach}$',
                   r'$tke/u_{\infty}^2*100$',
                   r'$\langle u^{\prime}u^{\prime}\rangle /u_{\infty}^2\cdot 100$',
                   r'$\langle u^{\prime}v^{\prime}\rangle /u_{\infty}^2\cdot 1000$',
                   r'$\sqrt{\langle p^{\prime}p^{\prime}\rangle} / p_{\infty}$']
    norm        = [params.u_ref, params.u_ref/100, params.u_ref/100, params.T_ref,
                   1.0, params.u_ref**2/100, params.u_ref**2/100, params.u_ref**2/1000,
                   params.p_ref]
    bbox        = [-100,20,-2,12,-20,20]
    streamline  = False
    D_norm      = False

# =============================================================================

    grid         = GridData( dirs.grid )
    grid.read_grid()
    blocklist, _ = grid.select_sliced_blockgrids('X', loc, bbox=bbox)

    stat_file    = os.path.join(dirs.sup_dir, stat_file)
    stat         = StatisticData(stat_file)
    stat.grid3d  = grid
    stat.read_statistic(blocklist,vars_in=vars_read)
    stat.match_grid(blocklist, grid, add_to_df=True )
    stat.compute_vars(blocklist,['mach','RS','p`'])
    
    roughwall = params.roughwall
    
    if roughwall:
        wdfile = get_filelist( dirs.wall_dist, 'snapshot.bin' )[0]
        wdsnap = Snapshot(wdfile)
        wdsnap.grid3d = grid
        wdsnap.read_snapshot(block_list=blocklist, var_read=['wd'])
        wdsnap = wdsnap.get_slice('X',loc, bbox=[-100,20,-2,0.5,0.0,5.2])
        wall   = wdsnap.get_contour('wd',0.0, blocklist=blocklist)
    else:
        wall   = np.array([[0,0,-100],[0,0,100]],dtype=float) 
    
    # periodic averaging

    stat.spanwise_periodic_average( blocklist, vars_output, params.D )

    dataset = pv.MultiBlock( stat.create_vtk_multiblock(blocklist,vars_output) )
    dataset = dataset.cell_data_to_point_data()

    # interpolate into a whole cartesian grid

    pz = np.linspace( 0.0, 5.2, 101,  endpoint=True)
    py = np.linspace(-0.6, 5.2, 201,  endpoint=True)
    px = np.array([0.0])
    df = pv_interpolate( dataset, vars_output, [px,py,pz] )
    
    os.chdir( create_folder(os.path.join(dirs.pp_statistics, outfolder)) )
    
    for i, var in enumerate(vars_output):
        output_var_fig( var, df, params, norm[i],  cbar_range[i], cbar_label[i],  wall=wall, D_norm=D_norm, streamline=streamline )
    

def output_var_fig(var, df, params:Params, norm, cbar_range, cbar_label, wall=None, D_norm=False, streamline=False):

    py = np.unique(df['y'])
    pz = np.unique(df['z'])

    data = np.array( df[var]   ).reshape( (len(py),len(pz)) )/norm
    mach = np.array( df['mach']).reshape( (len(py),len(pz)) )
    w    = np.array( df['w']   ).reshape( (len(py),len(pz)) )/params.u_ref
    v    = np.array( df['v']   ).reshape( (len(py),len(pz)) )/params.u_ref

    if D_norm: 
        length_unit = params.D
        y_bottom    = -0.48
        x_lim       = [0.0,1.0]
        y_lim       = [-0.48,1.6] 
        x_ticks     = [0.0,0.5,1.0]
        x_ticklabels= [r'$0.0$',r'$0.5$',r'$1.0$']
        y_ticks     = [-0.4,0.0,0.4,0.8,1.2,1.6]
        x_label     = r'$z/D$'
        y_label     = r'$y/D$'
        loc_tag     = [0.95, 1.4]
    else     : 
        length_unit = params.delta_0
        y_bottom    = -0.12
        x_lim       = [0.0,1.0]
        y_lim       = [-0.12,1.0]
        x_ticks     = [0.0,0.5,1.0]
        x_ticklabels= [r'$0.00$',r'$0.5$',r'$1.0$']
        y_ticks     = [-0.1,0.0,0.5,1.0]
        x_label     = r'$z/\delta_0$'
        y_label     = r'$y/\delta_0$'
        loc_tag     = [0.8, 0.8]
    
#    ywall = define_wall_shape(pz, casecode=params.casecode, write=False)/length_unit
    z     = np.array( df['z']    ).reshape( (len(py),len(pz)) )/length_unit
    y     = np.array( df['y']    ).reshape( (len(py),len(pz)) )/length_unit

    # visualization

    fig = plt.figure(figsize=(10,6))
    ax  = fig.add_axes([0.1,0.2,0.90,0.7])

    cbar_levels = np.linspace( cbar_range[0], cbar_range[1], 51)
    cbar_ticks  = np.linspace( cbar_range[0], cbar_range[1], 5)

    if var == 'u`v`': cmap = 'RdBu'
    else: cmap = 'RdBu_r'
    cs     = ax.contourf(z, y, data, levels=cbar_levels, cmap=cmap, extend='both')
    csnoic = ax.contour( z, y, mach, levels=[1.0], colors='lime', linewidths=2.0, zorder=9)

    if streamline:
        ax.streamplot(z,y, w, v, color='black', linewidth=0.5, density=1.0)
        #ax.quiver(z[::4,::4], y[::4,::4], w[::4,::4], v[::4,::4], color='black', scale=100)

#    ax.fill_between( pz/length_unit, y_bottom, y2=ywall, color='gray', zorder=10 )
    ax.fill_between( wall[:,2]/length_unit, y_bottom, y2=wall[:,1]/length_unit, color='gray', zorder=10 )
#    ax.plot( wall[:,2]/length_unit, wall[:,1]/length_unit, zorder=10 )
    ax.set_aspect('equal')

    cbar = plt.colorbar( cs, 
                        orientation='vertical', 
                        location='left', 
                        aspect=10,
                        ticks=cbar_ticks,
                        pad=0.30,
                        shrink=0.6) 

    cbar.outline.set_linewidth(1.5)

    cbar.ax.tick_params( direction='in',
                        left=True,right=False,
                        labelleft=True,labelright=False,
                        length=5,
                        width=1.0)

    cbar.ax.set_xlabel(cbar_label,labelpad=20)

    ax.set_xlim( x_lim )
    ax.set_ylim( y_lim )

    ax.set_xticks( x_ticks )
    ax.set_xticklabels( x_ticklabels )
    ax.set_yticks( y_ticks )
    ax.tick_params( which='major',
                    axis='both', 
                    direction='out',
                    length=10.0,
                    width=1.0)
    ax.tick_params(axis='y', pad=15)
    ax.tick_params(axis='x', pad=10)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.text( loc_tag[0], loc_tag[1],
            params.tag,
            va='center',
            ha='right',
            zorder=12,
            bbox={"fc":"white","alpha":0.8,"ec":"None"})    

    ax.spines[:].set_color('black')
    ax.spines[:].set_linewidth(1.5)
    ax.spines[:].set_zorder(11)

    figname     = 'zoomin_'+var
    if streamline: figname += '_streamline'
    if D_norm:     figname += '_D'
    figname    += '.png'

    plt.savefig(figname, dpi=300)
    # plt.show()
    plt.close()
    print(f"Figure saved to {os.getcwd()}/{figname}.")


if __name__ == "__main__":
    main()