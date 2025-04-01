#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_slicex.py
@Time    :   2023/09/11 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import time
import pickle
import numpy             as     np
from   scipy.interpolate import griddata

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.log         import Logger
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.tools       import define_wall_shape
from   vista.plane_analy import save_sonic_line
from   vista.plane_analy import shift_coordinates
from   vista.plane_analy import periodic_average
from   vista.plot_style  import plot_slicex_stat
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

locs_delta = np.linspace(-20,-20,1)
outfolder  = '/yz_planes'
periodic_ave = False

# =============================================================================

vars = ['u','v','w','uu','vv','ww','uv','vw','T','p','pp','rho']

dirs = Directories( os.getcwd() )

datafile = dirs.statistics
gridfile = dirs.grid
outpath  = dirs.pp_statistics + outfolder

# - read in case paramters

params  = Params( dirs.case_para_file )
delta    = params.delta_0
h_ridge  = params.H
h_md     = params.H_md
x_imp    = params.x_imp
p_ref    = params.p_ref
u_ref    = params.u_ref
casecode = params.casecode
n_period = params.n_period
tag      = params.tag

locs = locs_delta*delta + x_imp


# - read in grid info

G = GridData( gridfile )
G.read_grid()

# - enter outpath

create_folder( outpath )
os.chdir(outpath)


# - do slicing and output slices

for i, loc in enumerate(locs):

    loc_delta = locs_delta[i]
    df_slice_file = f"df_slice_{i:02d}.pkl"
    title = 'x= '+str(loc_delta)
    
    # check if the slice is already done
    
    if not os.path.exists(df_slice_file):
        
        print(f"Start doing slicing at x = {loc_delta:10.2f}.\n")

        block_list, indx_slic = G.select_sliced_blockgrids( 'X', loc )

        print(f"Selected {len(block_list)} blocks.\n")

        # - read statistics data file

        S = StatisticData( datafile )
        S.verbose = True
            
        S.read_statistic( block_list=block_list, vars_in=vars)
            
        S.match_grid( block_list, G )
        
        S.compute_vars( block_list, ['mach','RS','p`'])
        
        S.compute_gradients( block_list, ['vorticity','Q_cr','lambda2'])
        
        S.compute_source_terms( block_list, G )
                
        with timer("Get slice dataframe"):
            
            df_slice = S.get_slice_df( block_list, G, indx_slic, 'X' )
            
            with open(df_slice_file,'wb') as f:
                pickle.dump( df_slice, f )
    
    # - read in slice dataframe
    
    else: 
        print(f"{df_slice_file} already exists, read in directly...\n")
        df_slice = pickle.load( open(df_slice_file,'rb') )
    
    
# ----- interpolate and plot    
        
    with timer("Interpolate and plot "):
               
        df_slice = shift_coordinates( df_slice, delta, h_ridge, h_md, x_imp )  
        
        y_slice = np.array( df_slice['ys'] )
        z_slice = np.array( df_slice['zs'] )
        
        mach_slice = np.array( df_slice['mach'] )
        tke_slice = np.array( df_slice['tke'] )
        RS_slice = np.array( df_slice['u`v`'] )
        uu_slice = np.array( df_slice['u`u`'] )
        S_slice = np.array( df_slice['S'] )
        w1_slice = np.array( df_slice['w1'] )
        u_slice = np.array( df_slice['u'] )
        v_slice = np.array( df_slice['v'] )
        w_slice = np.array( df_slice['w'] )
        p_fluc_slice = np.array( df_slice['p`'] )
        rho_slice = np.array( df_slice['rho'] )
        Q_cr_slice = np.array( df_slice['Q_cr'] )
        l2_slice   = np.array( df_slice['lambda2'] )
        ke_slice   = 0.5*(v_slice**2 + w_slice**2)
        
        # generate interpolation grid
        
        z = np.linspace(-1.0,1.0, 320)  # must be muliple of 16 for periodic average
        if 'smooth' in casecode:
            y = np.linspace(0.02, 1.1, 55)
        else:
            y = np.linspace(-0.1, 1.1, 241)
        
        zz,yy = np.meshgrid(z,y)
        
        # mapping variables
        
        mach = griddata( (z_slice,y_slice), mach_slice,
                         (zz,yy), method='linear')
        
        tke = griddata( (z_slice,y_slice), tke_slice,
                        (zz,yy), method='linear')

        RS = griddata( (z_slice,y_slice), RS_slice,
                        (zz,yy), method='linear')
        
        uu = griddata( (z_slice,y_slice), uu_slice,
                        (zz,yy), method='linear')
        
        S = griddata( (z_slice,y_slice), S_slice,
                      (zz,yy), method='linear')
        
        w1 = griddata( (z_slice,y_slice), w1_slice,
                      (zz,yy), method='linear')    

        u = griddata( (z_slice,y_slice), u_slice,
                      (zz,yy), method='linear') 

        v = griddata( (z_slice,y_slice), v_slice,
                      (zz,yy), method='linear') 

        w = griddata( (z_slice,y_slice), w_slice,
                      (zz,yy), method='linear')

        p_fluc = griddata( (z_slice,y_slice), p_fluc_slice,
                           (zz,yy), method='linear')   

        rho = griddata( (z_slice,y_slice), rho_slice,
                        (zz,yy), method='linear')  
        
        Q_cr = griddata( (z_slice,y_slice), Q_cr_slice,
                         (zz,yy), method='linear')
    
        l2 = griddata( (z_slice,y_slice), l2_slice,
                       (zz,yy), method='linear')
        
        ke = griddata( (z_slice,y_slice), ke_slice,
                       (zz,yy), method='linear')
        
# ------ extending corner grid for smooth wall

        if 'smooth' in casecode:
            
            len_ext = np.shape(zz)[1]
            zz = np.concatenate(([zz[0,:]],zz),axis=0)
            yy = np.concatenate(([np.zeros(len_ext)],yy),axis=0)
            mach = np.concatenate(([np.zeros(len_ext)],mach),axis=0)
            tke = np.concatenate(([np.zeros(len_ext)],tke),axis=0)
            RS = np.concatenate(([np.zeros(len_ext)],RS),axis=0)
            uu = np.concatenate(([np.zeros(len_ext)],uu),axis=0)
            S  = np.concatenate(([np.zeros(len_ext)],S),axis=0)
            w1 = np.concatenate(([w1[0,:]],w1),axis=0)
            u  = np.concatenate(([np.zeros(len_ext)],u),axis=0)
            v  = np.concatenate(([np.zeros(len_ext)],v),axis=0)
            w  = np.concatenate(([np.zeros(len_ext)],w),axis=0)
            p_fluc = np.concatenate(([p_fluc[0,:]],p_fluc),axis=0)
            Q_cr = np.concatenate(([np.zeros(len_ext)],Q_cr),axis=0)
            l2 = np.concatenate(([np.zeros(len_ext)],l2),axis=0)
            ke = np.concatenate(([np.zeros(len_ext)],ke),axis=0)
        
        
        save_sonic_line( zz, yy, mach )
        
        define_wall_shape( z*5.2, casecode=casecode, yshift=(h_ridge-h_md) )
        
        # Do periodic average for smallest ridge spacing case
        
        if periodic_ave:  # otherwise streamline looks too messy
            n_period = int(n_period/2)
            w = periodic_average(w,n_period,axis=1)
            v = periodic_average(v,n_period,axis=1)
            mach   = periodic_average(mach,n_period,axis=1)
            tke    = periodic_average(tke,n_period,axis=1)
            RS     = periodic_average(RS,n_period,axis=1)
            uu     = periodic_average(uu,n_period,axis=1)
            w1     = periodic_average(w1,n_period,axis=1)
            p_fluc = periodic_average(p_fluc,n_period,axis=1)
            rho    = periodic_average(rho,n_period,axis=1)
            Q_cr   = periodic_average(Q_cr,n_period,axis=1)
            l2     = periodic_average(l2,n_period,axis=1)
            ke     = periodic_average(ke,n_period,axis=1)

# ------ Mach number and streamlines

        cbar = r'$\langle Mach \rangle$'
        cbar_levels = np.linspace(0.0,2.0,21)
        cbar_ticks  = np.linspace(0.0,2.0,5)
        if 'smooth' in casecode:
            plot_slicex_stat( zz, yy, mach,
                              tag=tag,
                              filename=casecode+'_Mach_'+str(i+1),
                              col_map='RdBu_r',
                              cbar_label=cbar,
                              cbar_levels=cbar_levels,
                              cbar_ticks=cbar_ticks,
                              title=title)
        else:
            plot_slicex_stat( zz, yy, mach,
                              vectors=[rho*w,rho*v],
                              tag=tag,
                              filename=casecode+'_Mach_'+str(i+1),
                              col_map='RdBu_r',
                              cbar_label=cbar,
                              cbar_levels=cbar_levels,
                              cbar_ticks=cbar_ticks,
                              title=title)

# ------ source term

#        cbar = r'$S$'
#        cbar_levels=np.linspace(-20000.0,20000.,51)
#        plot_slicex_stat( zz, yy, S,
#                          tag=tag,
#                          filename='S_'+str(i+1),
#                          cbar_label=cbar,
#                          cbar_levels=cbar_levels,
#                          title=title)

# ------ vorticity
  
        cbar = r'$\frac{\langle \omega_x \rangle \delta_0}{u_\infty}$ '
        cbar_levels = np.linspace(-1.0,1.0,26)
        cbar_ticks  = np.linspace(-1.0,1.0,5)
        plot_slicex_stat( zz, yy, w1*delta/u_ref,
                          tag=tag,
                          vectors=[w,v],
                          arrow=True,
                          filename=casecode+'_vorticity_'+str(i+1),
                          col_map='RdBu_r',
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          title=title)

# ------ streamwise velocity u

        cbar = r'$\langle u \rangle /u_{\infty}$'
        cbar_levels = np.linspace(0.0,1.0,21)
        cbar_ticks  = np.linspace(0.0,1.0,5)
        plot_slicex_stat( zz, yy, u/u_ref,
                          tag=tag,
                          filename=casecode+'_u_'+str(i+1),
                          cbar_label=cbar,
                          col_map='RdBu_r',
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          title=title)

# ------ vertical velocity v

        cbar = r'$\frac{\langle v \rangle  \cdot 100}{u_{\infty}}$'
        cbar_levels = np.linspace(-3.0,3.0,26)
        cbar_ticks  = np.linspace(-3.0,3.0,5)
        plot_slicex_stat( zz, yy, v/u_ref*100,
                          tag=tag,
                          vectors=[w,v],
                          arrow=True,
                          filename=casecode+'_v_'+str(i+1),
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          col_map='RdBu_r',
                          title=title)

        print(f"vmax = {np.max(v):10.2f} (vmax/u_ref*100 ={np.max(v/u_ref*100):10.2f})")

# ------ spanwise velocity w

        cbar = r'$\frac{\langle w \rangle  \cdot 100}{u_{\infty}}$'
        cbar_levels = np.linspace(-2.0,2.0,26)
        cbar_ticks  = np.linspace(-2.0,2.0,5)
        plot_slicex_stat( zz, yy, w/u_ref*100,
                          tag=tag,
                          vectors=[w,v],
                          arrow=True,
                          filename=casecode+'_w_'+str(i+1),
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          col_map='RdBu_r',
                          title=title)

        print(f"wmax = {np.max(w):10.2f} (wmax/u_ref*100 ={np.max(w/u_ref*100):10.2f})")


# ------ turbulent kinetic energy

        cbar = r'$\langle tke ^+ \rangle \cdot 100$'
        cbar_levels = np.linspace(0.0,2.5,26)
        cbar_ticks  = np.linspace(0.0,2.5,6)
        plot_slicex_stat( zz, yy, tke/(u_ref**2)*100,
                          tag=tag,
                          filename=casecode+'_tke_'+str(i+1),
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          col_map='RdBu_r',
                          title=title)

# ------ Reynolds stresses

        cbar = r"$-\langle u^{'} v^{'} \rangle$"
        plot_slicex_stat( zz, yy, -RS/(u_ref**2),
                          tag=tag,
                          filename=casecode+'_RS_'+str(i+1),
                          cbar_label=cbar,
                          col_map='RdBu_r',
                          title=title)

        cbar = r"$-\langle u^{'} u^{'} \rangle$"
        cbar_levels = np.linspace(0.0,2.0,26)
        cbar_ticks  = np.linspace(0.0,2.0,6)
        plot_slicex_stat( zz, yy, uu/(u_ref**2)*100,
                          tag=tag,
                          filename=casecode+'_uu_'+str(i+1),
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          col_map='RdBu_r',
                          title=title)
        
# ----- pressure fluctuation

        cbar = r"$\sqrt{\langle p^{'}p^{'}\rangle}\cdot 100$"
        cbar_levels = np.linspace(0.0, 3.2,33)
        cbar_ticks  = np.linspace(0.0,3.2,5)
        plot_slicex_stat( zz, yy, p_fluc/p_ref*100,
                          tag=tag,
                          filename=casecode+'_p_fluc_'+str(i+1),
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          col_map='RdBu_r',
                          cbar_ticks=cbar_ticks,
                          title=title)

# ------ Q_cr
  
        cbar = r'$Q$'
        Q_max = np.max(Q_cr)
        cbar_levels = np.linspace(0,800,25)
#        cbar_ticks  = np.linspace(-1.0,1.0,5)
        plot_slicex_stat( zz, yy, Q_cr,
                          tag=tag,
                          vectors=[w,v],
                          arrow=True,
                          filename=casecode+'_Q_'+str(i+1),
                          col_map='RdBu_r',
                          cbar_label=cbar,
                          wall=False,
                          cbar_levels=cbar_levels,
#                          cbar_ticks=cbar_ticks,
                          title=title)

# ------ lambda2
  
        cbar = r'$Q$'
        l2_max = np.max(l2)
        cbar_levels = np.linspace(-500,0,25)
#        cbar_ticks  = np.linspace(-1.0,1.0,5)
        plot_slicex_stat( zz, yy, l2,
                          tag=tag,
                          vectors=[w,v],
                          arrow=True,
                          filename=casecode+'_l2_'+str(i+1),
                          col_map='RdBu_r',
                          cbar_label=cbar,
                          wall=False,
                          cbar_levels=cbar_levels,
#                          cbar_ticks=cbar_ticks,
                          title=title)    
        
# ------ kinetic energy
  
        cbar = r'$\langle k \rangle$'
        ke_max = np.max(ke)
        cbar_levels = np.linspace(0,0.5,26)
        
        plot_slicex_stat( zz, yy, ke,
                          tag=tag,
                          vectors=[w,v],
                          arrow=True,
                          filename=casecode+'_ke_'+str(i+1),
                          cbar_label=cbar,
                    #      cbar_levels=cbar_levels,
                    #      cbar_ticks=cbar_ticks,
                          col_map='plasma',
                          title=title)
        

    print(f"Finished doing slicing at x = {loc_delta:10.2f} delta.",end='') 
    print(f" Progress {(i+1)/len(locs)*100:5.2f} %.")
    print(f"=====================================\n")

# print out the time finishing the job
    
print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")