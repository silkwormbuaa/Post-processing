#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   show_slice.py
@Time    :   2023/09/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   get a slice from statistics.bin (only applicable to 3D)
             df_slice_00 is at the ridge top, df_slice_01 is at the valley.
             Be careful with the interpolation range for different ridge height 
             cases!!!
'''

import os
import sys
import time
import pickle
import numpy             as     np
from   scipy.interpolate import griddata

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.plane_analy import save_sonic_line
from   vista.plane_analy import save_isolines
from   vista.plane_analy import shift_coordinates
from   vista.plane_analy import compute_DS
from   vista.plot_style  import plot_slicez_stat
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

# locs = [ridge location, valley location]

outfolder = '/xy_planes'
vars      = ['u','v','w','T','rho','uu','vv','ww','uv','pp','p']

# =============================================================================

dirs = Directories( os.getcwd() )

datafile = dirs.statistics
gridfile = dirs.grid
outpath  = dirs.pp_statistics + outfolder

# - read in case params

params = Params( dirs.case_para_file )
delta   = params.delta_0
h_ridge = params.H
h_md    = params.H_md
x_imp   = params.x_imp
p_ref   = params.p_ref
u_ref   = params.u_ref
x_pfmax = params.x_pfmax
casecode =  params.casecode
target_dir = '/home/wencan/temp/DataPost/contour/'+casecode+'/'
create_folder( target_dir )

locs = [ 0.001, 0.5*params.D + 0.001 ]

# - read in grid info

G = GridData( gridfile )

G.read_grid()

# - enter outpath

create_folder( outpath )
os.chdir(outpath)

# - do slicing and output slices

for i, loc in enumerate( locs ):
    
    # check if the slice is already there
    
    df_slice_file = f"df_slice_{i:02d}.pkl"
    
    if os.path.exists( df_slice_file ):
        print(f"File {df_slice_file} already exists, read in directly...\n")
        df_slice = pickle.load( open( df_slice_file, 'rb' ) )
    
    else: # if not, do slicing
        block_list, indx_slic = G.select_sliced_blockgrids( 'Z', loc )

        print(f"Selected {len(block_list)} blocks.\n")

        # - read statistics data file

        S = StatisticData( datafile )
        S.grid3d = G

        with timer("read selected blocks "):
            
            S.read_statistic( block_list, vars )
            S.match_grid( block_list, G )
            S.compute_vars( block_list, ['mach','RS','p`'] )
            S.compute_gradients( block_list, ['grad_rho','shadowgraph','grad_p'])
            
        with timer("Get slice dataframe "):
            
            df_slice = S.get_slice_df( block_list, G, indx_slic, 'Z' )
            
            df_slice = shift_coordinates( df_slice, delta, h_ridge, h_md, x_imp )  
            
            with open( df_slice_file, 'wb' ) as f:
                pickle.dump( df_slice, f )
                
            os.system(f"cp {df_slice_file} {target_dir}")

    with timer("Interpolate and plot "):
        
# ----- interpolate
        
        x_slice = np.array( df_slice['xs'] )
        y_slice = np.array( df_slice['y_scale'] )
        
        u_slice       = np.array( df_slice['u'] )
        mach_slice    = np.array( df_slice['mach'] )
        gradrho_slice = np.array( df_slice['grad_rho'] )
        T_slice       = np.array( df_slice['T'] )
        tke_slice     = np.array( df_slice['tke'] )
        p_fluc_slice  = np.array( df_slice['p`'] )
        grad_p_slice  = np.array( df_slice['grad_p'] )
        
        x = np.linspace(-20,10,301)
        
        if i == 0:
            y = np.linspace(0.005, 8, 201)
            wall = False
        
        if i == 1:
            if casecode == 'smooth':  y = np.linspace(0.005, 8, 201)
            else:                     y = np.linspace(-0.10,8, 405)
            wall=False
        
        xx,yy = np.meshgrid(x,y)
        
        mach    = griddata( (x_slice,y_slice), mach_slice,
                            (xx,yy), method='linear')
        
        u       = griddata( (x_slice,y_slice), u_slice,
                            (xx,yy), method='linear')
        
        gradrho = griddata( (x_slice,y_slice), gradrho_slice,
                            (xx,yy), method='linear')

        T       = griddata( (x_slice,y_slice), T_slice,
                            (xx,yy), method='linear')

        tke     = griddata( (x_slice,y_slice), tke_slice,
                            (xx,yy), method='linear')

        p_fluc  = griddata( (x_slice,y_slice), p_fluc_slice,
                            (xx,yy), method='linear') 

        grad_p  = griddata( (x_slice,y_slice), grad_p_slice,
                            (xx,yy), method='linear') 


# ----- save sonic line, separation line, shock shape, DS

        soniclinefile = f'soniclines_{i:02d}.pkl'
        save_sonic_line( xx,yy, mach, out_file=soniclinefile )
        os.system(f"cp {soniclinefile} {target_dir}")
        
        seplinefile = f'seplines_{i:02d}.pkl'
        save_isolines( xx, yy, u, 0.0, seplinefile )
        os.system(f"cp {seplinefile} {target_dir}")
        
        shockshapefile = f'shockshape_{i:02d}.pkl'
        save_isolines( xx, yy, gradrho, 0.15, shockshapefile, clip=True)
        os.system(f"cp {shockshapefile} {target_dir}")
        
        shockDSfile = f'shockDS_{i:02d}.pkl'
        DS = compute_DS( gradrho, min=0.0, max=2.0 )
        save_isolines( xx, yy, DS, 0.2, shockDSfile, clip=True)
        os.system(f"cp {shockDSfile} {target_dir}")
        
        print(f"Saved {soniclinefile}, {seplinefile}, {shockshapefile}, {shockDSfile}.\n")

# ----- plot slices

        cbar = r'$Mach$'
        
        plot_slicez_stat( xx,yy,mach, 
                          separation=seplinefile,
                          shockshape=shockshapefile,
                          filename=casecode+'_MachZ_'+str(i+1),
                          cbar_label=cbar,
                          col_map='coolwarm',
                          wall=wall,
                          tag = params.tag,
                          tag_loc=[7.0,4.5],
                          x_lim=[-15,10],
                          y_lim=[0,6] )

        cbar = r'$<T>/T_{\infty}$'
        
        plot_slicez_stat( xx,yy,T/160.15,
                          separation=seplinefile,
                          shockshape=shockshapefile,
                          sonic=soniclinefile,
                          filename=casecode+'_TemperatureZ_'+str(i+1),
                          col_map='plasma',
                          cbar_label=cbar,
                          wall=wall,
                          tag = params.tag,
                          tag_loc=[7.0,4.5],
                          x_lim=[-15,10],
                          y_lim=[0,6] )

        cbar = r'$\nabla{\rho}$'
        
        plot_slicez_stat( xx,yy,gradrho, 
                          separation=seplinefile,
                          shockshape=shockshapefile,
                          sonic=soniclinefile,
                          filename=casecode+'_grad_rho_'+str(i+1),
                          col_map='Greys',
                          cbar_label=cbar,
                          wall=wall,
                          tag = params.tag,
                          tag_loc=[7.0,4.5],
                          x_lim=[-15,10],
                          y_lim=[0,6] )

        plot_slicez_stat( xx,yy,gradrho, 
                          separation=seplinefile,
                          shockshape=shockshapefile,
                          sonic=soniclinefile,
                          filename=casecode+'_mean_shock_shape_'+str(i+1),
                          col_map='Greys',
                          cbar_label=cbar,
                          wall=wall,
                          tag = params.tag,
                          tag_loc=[7.0,4.5],
                          x_lim=[-15,10],
                          y_lim=[0,6] )

        cbar = 'DS'
        plot_slicez_stat( xx,yy,DS, 
                          separation=seplinefile,
                          DS=shockDSfile,
                          sonic=soniclinefile,
                          filename=casecode+'_mean_shock_DS_'+str(i+1),
                          col_map='Greys_r',
                          cbar_label=cbar,
                          wall=wall,
                          tag = params.tag,
                          tag_loc=[7.0,4.5],
                          x_lim=[-15,10],
                          y_lim=[0,6] )


        cbar = r'$tke$'
        cbar_levels = np.linspace(0,12500,51)
        cbar_ticks  = np.linspace(0,12500,5)
        
        plot_slicez_stat( xx,yy,tke, 
                          separation=seplinefile,
                          shockshape=shockshapefile,
                          sonic=soniclinefile,
                          x_pfmax=x_pfmax,
                          filename=casecode+'_tke_'+str(i+1),
                          col_map='coolwarm',
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          wall=wall,
                          tag = params.tag,
                          tag_loc=[7.0,4.5],
                          x_lim=[-15,10],
                          y_lim=[0,6] )
        
        cbar = r"$\sqrt{\langle p'p' \rangle}/p_{\infty}$"
        cbar_levels = np.linspace(0,0.5,51)
        cbar_ticks  = np.linspace(0,0.5,6)
        
        plot_slicez_stat( xx,yy,p_fluc/p_ref, 
                          separation=seplinefile,
                          shockshape=shockshapefile,
                          sonic=soniclinefile,
                          x_pfmax=x_pfmax,
                          filename=casecode+'_pressure_fluc_'+str(i+1),
                          col_map='coolwarm',
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          wall=wall,
                          tag = params.tag,
                          tag_loc=[7.0,4.5],
                          x_lim=[-15,10],
                          y_lim=[0,6] )

        cbar = 'pressure gradient'
        cbar_levels = np.linspace(0,12,51)
        cbar_ticks  = np.linspace(0,12,7)
        plot_slicez_stat( xx,yy,grad_p*delta/p_ref, 
                          separation=seplinefile,
                          shockshape=shockshapefile,
                          sonic=soniclinefile,
                          x_pfmax=x_pfmax,
                          filename=casecode+'_p_grad_'+str(i+1),
                          col_map='coolwarm',
                          cbar_label=cbar,
                          cbar_levels=cbar_levels,
                          cbar_ticks=cbar_ticks,
                          wall=wall,
                          tag = params.tag,
                          tag_loc=[7.0,4.5],
                          x_lim=[-15,10],
                          y_lim=[0,6] )
        
    
# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    
sys.stdout.flush()       
    

