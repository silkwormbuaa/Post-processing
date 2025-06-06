#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   spanwise_average.py
@Time    :   2025/03/14 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   just do direct ensemble average for spanwise direction,
             (no fluid/solid cells are considered).
             Results are save into a pickled pyvista dataset.
'''

import os
import sys
import numpy              as     np
import pyvista            as     pv
import matplotlib.pyplot  as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid         import GridData
from   vista.statistic    import StatisticData
from   vista.directories  import Directories
from   vista.params       import Params
from   vista.tools        import lin_grow
from   vista.directories  import create_folder
from   vista.plane_analy  import save_isolines
from   vista.plane_analy  import pv_interpolate
from   vista.plot_setting import set_plt_rcparams

def main():

    case_dir   = '/home/wencan/temp/smooth_adiabatic/'
    vars_read  = ['u','v','w','p','T','pp','uu','vv','ww','uv']
    vars_out   = ['u','v','w','p','T','mach','grad_rho','tke','u`u`','v`v`','w`w`','u`v`','p`','Px','Py','Ps','Pk']
    vars_load  = ['u','mach','grad_rho','tke','u`u`','v`v`','w`w`','u`v`','p`','Px','Py','Ps','Pk']    
    bbox       = [-57.5,120,0,100,0.0,20]
    rescale    = [-50.4, 0.0, 0.0, 5.2, 5.2, 5.2]
    clean      = False
    
    dirs       = Directories( case_dir )
    params     = Params( dirs.case_para_file )
    
    os.chdir( create_folder(dirs.pp_z_average) )

    vars_show  = ['u`u`','v`v`','w`w`','tke','p`','Px','Py','Ps','Pk']
    cbar_lvls  = [[0,9_000],[0,3_000],[0,5_000],[0,15_000],[0,0.15],[0,400_000],[0,400_000],[0,400_000],[0,400_000]]
    labels     = [r"$u'u'$",r"$v'v'$",r"$w'w'$",'tke',r"$\sqrt{\langle p'p' \rangle}/p_{\infty}$",r'$P_x$',r'$P_y$',r'$P_s$',r'$P_k$']
    norms      = [1.0,1.0,1.0,1.0,params.p_ref,1.0,1.0,1.0,1.0,1.0]
    
    dataset = load_dataset( dirs, bbox, vars_read, vars_out, rescale, clean )
    set_plt_rcparams(latex=True,fontsize=25)
    
    os.chdir( create_folder(dirs.pp_z_average) )
    for i, var in enumerate(vars_show):
        post_process_dataset( params, dataset, vars_load, var, cbar_lvls[i], norms[i], labels[i] )
#    visualize( dataset, 'u`u`' )

# =============================================================================

def load_dataset( dirs:Directories, bbox, vars_read, vars_out, rescale, clean=False):
    
    if clean: os.system( f'rm -rf {dirs.pp_z_average}/z_average*' )
    
    if not os.path.exists( dirs.pp_z_average + '/z_average.vtmb' ):

        print("z_average.vtmb not found, start to compute spanwise average...\n")
        
        create_folder( dirs.pp_z_average )
        
        grid = GridData( dirs.grid ) 
        grid.read_grid()
        blocklist = grid.select_blockgrids( bbox=bbox, mode='within' )

        stat = StatisticData( dirs.statistics )
        stat.grid3d = grid
        stat.read_statistic(    blocklist, vars_read )
        stat.compute_vars(      blocklist, vars_new=['mach','p`'] )
        stat.compute_gradients( blocklist, grads=['grad_rho'] )
        stat.compute_vars(      blocklist, vars_new=['RS','production'])
        dataset = pv.MultiBlock(stat.spanwise_average( blocklist, vars_out, rescale=rescale ))
        dataset.save( dirs.pp_z_average + '/z_average.vtmb', binary=True )

    else:
        
        print("Found z_average.vtmb, loading it...\n")
        
        dataset = pv.read( dirs.pp_z_average + '/z_average.vtmb' )
    
    return dataset

# =============================================================================

def post_process_dataset( params:Params, dataset:pv.MultiBlock, vars_load, var, lvls=None, norm=1.0, label=None):

    px   = np.linspace( -15, 10, 301, endpoint=True )
    py,_ = lin_grow( 0.0, 0.01, 1.04, upbound=11.0 )
    
    # points near interface where mesh density jumps will lead to holes
    
    # px   = np.linspace( -4.2, -3.4, 11,  endpoint=True ) 
    # py   = np.linspace(  3.3,  3.5, 101, endpoint=True )
    pz   = np.array([0.0])

    dataset = dataset.cell_data_to_point_data().combine()

    df = pv_interpolate( dataset, vars_load,[px,py,pz])

    x        = np.array( df['x']    ).reshape( (len(py),len(px)) )
    y        = np.array( df['y']    ).reshape( (len(py),len(px)) )
    u        = np.array( df['u']    ).reshape( (len(py),len(px)) )
    mach     = np.array( df['mach'] ).reshape( (len(py),len(px)) )
    grad_rho = np.array( df['grad_rho'] ).reshape( (len(py),len(px)) )
    var_data = np.array( df[var] ).reshape( (len(py),len(px)) ) / norm

    save_isolines( x, y, u,        0.0,  'sepline.pkl'   )
    save_isolines( x, y, mach,     1.0,  'sonicline.pkl' )
    save_isolines( x, y, grad_rho, 0.15, 'shockshape.pkl')

    fig, ax = plt.subplots(1,1,figsize=(12.8,7.2))
    if lvls is None:
        c    = ax.contourf( x, y, var_data, cmap='coolwarm', extend='both')
    else:
        c    = ax.contourf( x, y, var_data, levels=np.linspace(lvls[0],lvls[1],51), cmap='coolwarm', extend='both')
    csep = ax.contour(  x, y, u,        levels=[0.0],  colors='red',   linewidths=3.0, zorder=10)
    cson = ax.contour(  x, y, mach,     levels=[1.0],  colors='lime',  linewidths=3.0, zorder=10)
    cshk = ax.contour(  x[20:,:], y[20:,:], grad_rho[20:,:], 
                        levels=[0.15], colors='black', linewidths=2.0, zorder=9)
    ax.plot( params.x_pfmax, 0.0, '*', color='black', ms=20, zorder=11)
    
    ax.set_aspect('equal')
    
    ax.set_xlabel(r'$(x-x_{imp})/\delta_0$')
    ax.set_ylabel(r'$y/\delta_0$')
    ax.text(-12.5, 2.5, params.tag, fontsize=30 )
    
    ax.spines[:].set_linewidth(2.0)
    ax.tick_params( direction='out', length=10, width=2.0)
    ax.set_xlim([-13,-2.5])
    ax.set_ylim([0,3])
    
    cbar = fig.colorbar( c, 
                         ax=ax, 
                         pad=0.20,
                         shrink=0.5,
                         orientation='horizontal')
    if lvls is not None:
        cbar.set_ticks( np.linspace(lvls[0],lvls[1],6) )
        
    cbar.ax.tick_params( direction='in',
                         length=5,
                         width=1.5)
    cbar.ax.set_xlabel(label, loc='center', labelpad=20)
    cbar.outline.set_linewidth(1.5)
    
    plt.savefig(var+'.png', dpi=300, bbox_inches='tight')

    plt.close()
    
    print(f"Output {var}.png to {os.getcwd()}")

# =============================================================================

def visualize( dataset, var ):
    p = pv.Plotter()
    p.add_mesh( dataset, scalars=var, show_edges=False, cmap='coolwarm')
    p.add_axes()
    p.view_vector([0.0,0.0,1.0], viewup=[0.0,1.0,0.0])
    p.show()
    p.close()


# =============================================================================

if __name__ == '__main__':
    main()

