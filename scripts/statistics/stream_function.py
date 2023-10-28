#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   stream_function.py
@Time    :   2023/10/27 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import time
import pickle

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

import numpy             as     np

import pandas            as     pd

import matplotlib.pyplot as     plt

from   scipy.interpolate import griddata

from   vista.statistic   import StatisticData

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.tools       import define_wall_shape
from   vista.tools       import read_case_parameter
from   vista.tools       import get_filelist
from   vista.tools       import lin_grow

from   vista.plane_analy import shift_coordinates
from   vista.plane_analy import periodic_average
from   vista.plane_analy import compute_stream_function

from   vista.plot_style  import plot_slicex_stat

from   vista.log         import Logger
sys.stdout = Logger()
# =============================================================================

loc_delta = -20.0
outfolder  = '/stream_function'
periodic_ave = True

# =============================================================================

datapath = os.getcwd()

datafile = datapath + '/statistics.bin'
gridfile = datapath + '/inca_grid.bin'
ccfile   = datapath + '/cutcells_setup.dat'
outpath  = datapath + outfolder
parametersfile = datapath.split('/results')[0] + '/case_parameters'

# - read in case paramters

parameters = read_case_parameter( parametersfile )
delta   = float( parameters.get('delta_0') )
h_ridge = float( parameters.get('H') )
h_md    = float( parameters.get('H_md') )
x_imp   = float( parameters.get('x_imp') )
p_ref   = float( parameters.get('p_ref') )
u_ref   = float( parameters.get('u_ref') )
casecode =  str( parameters.get('casecode') )
n_period = int( parameters.get('period') )
tag      = str( parameters.get('tag'))
roughwall  = True if parameters.get('roughwall').lower() == 'true' else False

if roughwall:
    snapshotfile = get_filelist(datapath.split('/results')[0] +'/wall_dist',
                                key='snapshot.bin')[0]

loc = loc_delta*delta + 50.4


# - enter outpath

if not os.path.exists(outpath): 
    os.mkdir( outpath )
    print(f"Created directory {outpath}.\n")

os.chdir(outpath)

# - read in grid info

G = GridData( gridfile )

G.read_grid()

block_list, indx_slic = G.select_sliced_blockgrids( 'X', loc )

# - assign vol_fra to grid, then match G to S

with timer("Assign vol_fra to G"):
    
    if roughwall:

        # - read in wall distance data
        with timer("read wall distance field"):
            
            wd_snap = Snapshot( snapshotfile )
            wd_snap.read_snapshot( block_list )

        # - read in cut cell data and assign vol_fra

        with timer("read in cut cell info "):

            cc_df = pd.read_csv( ccfile, delimiter = r'\s+' )

            cc_df.drop( columns=['fax0','faz0'
                                ,'fax1','faz1', 'processor']
                                , inplace=True )
            
            for num in block_list:

                # dataframe slice for a certain block
                temp_df = cc_df[ cc_df['block_number'] == num ]
                
                wall_dist = np.array( wd_snap.snap_data[num-1].df['wd'] )
                
                # block number starts from 1, but python list index
                # starts from 0
                G.g[num-1].assign_vol_fra( df=temp_df, wall_dist=wall_dist )
    else:
        
        for num in block_list:
            G.g[num-1].assign_vol_fra()


# - do slicing and output slices

df_slice_file = f"df_slice_streamfunction.pkl"
title = 'x= '+str(loc_delta)

# check if the slice is already done

if not os.path.exists(df_slice_file):
    
    print(f"Start doing slicing at x = {loc_delta:10.2f}.\n")

    print(f"Selected {len(block_list)} blocks.\n")
    
    # read statistics data file

    S = StatisticData( datafile )

    with timer("read selected blocks and match grids "):
        
        with open(datafile,'br') as f:
            
            S.read_stat_header( f )
            
            vars = ['rho','urho','vrho','wrho']
            
            S.read_stat_body( f, block_list, vars )
            
        S.compute_vars( block_list, ['favre_velocity'])
    
        # match grid and pass vol_fra from G to S
        
        S.match_grid( G, block_list )
        
        # manipulate rho (rho*vol_fra)
        
        for num in block_list:
            
            S.bl[num-1].df['rho'] = (np.array( S.bl[num-1].df['rho'] ) * 
                                     np.array( S.bl[num-1].df['vol_fra'] ) )
            
            i_belowwall = np.where(np.array( S.bl[num-1].df['vol_fra'])<1.0)[0]
            rho_slice = np.array( S.bl[num-1].df['rho'] )
            w_favre_slice = np.array( S.bl[num-1].df['w_favre'] )
            v_favre_slice = np.array( S.bl[num-1].df['v_favre'] )
            rho_slice[i_belowwall] = 0.0
            w_favre_slice[i_belowwall] = 0.0
            v_favre_slice[i_belowwall] = 0.0           
            S.bl[num-1].df['rho'] = rho_slice
            S.bl[num-1].df['w_favre'] = w_favre_slice
            S.bl[num-1].df['v_favre'] = v_favre_slice
            
    with timer("Get slice dataframe and match grids"):
        
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
    
    rho_slice  = np.array( df_slice['rho'] )
    w_favre_slice = np.array( df_slice['w_favre'] )
    v_favre_slice = np.array( df_slice['v_favre'] )
    
    # generate interpolation grid
    
    z = np.linspace(-1.0,1.0, 320)
    if casecode == 'smooth_wall':
        y = np.linspace(0.02, 1.1, 55)
    else:
        y = np.linspace(-0.1, 1.1, 481)
    
    zz,yy = np.meshgrid(z,y)
    
    # mapping variables
    
    rho = griddata( (z_slice,y_slice), rho_slice,
                    (zz,yy), method='linear')
    
    w_favre = griddata( (z_slice,y_slice), w_favre_slice,
                        (zz,yy), method='linear')
    
    v_favre = griddata( (z_slice,y_slice), v_favre_slice,
                        (zz,yy), method='linear')

    # Do periodic average for smallest ridge spacing case
    
    if periodic_ave:  # otherwise streamline looks too messy
        n_period = int(n_period/2)
        rho = periodic_average(rho,n_period,axis=1,sym=True)  
        w_favre = periodic_average(w_favre,n_period,axis=1,antisym=True)  
        v_favre = periodic_average(v_favre,n_period,axis=1,sym=True)    
    
# -- compute stream function
    
    rho_w = 0.58050
    
    psi = compute_stream_function( rho, rho_w, w_favre, v_favre, z, y )


# -- extending corner grid for smooth wall

    if casecode == 'smooth_wall':
        
        len_ext = np.shape(zz)[1]
        zz = np.concatenate(([zz[0,:]],zz),axis=0)
        yy = np.concatenate(([np.zeros(len_ext)],yy),axis=0)
        psi = np.concatenate(([np.zeros(len_ext)],psi),axis=0)
    
    
# -- find the extreme value of stream function and identify the location

    psi_max = np.max(psi)   
    jmax,kmax = np.where(psi==psi_max)
    y_max = yy[jmax,kmax]; z_max = zz[jmax,kmax]
    
    psi_min = np.min(psi)   
    jmin,kmin = np.where(psi==psi_min)
    y_min = yy[jmin,kmin]; z_min = zz[jmin,kmin]
    
    print(f"Stream function maximum {psi_max:10.5f} at ",end='')
    print(f"z = {z_max[0]:10.5f} and y = {y_max[0]:10.5f}.\n")
    print(f"Stream function minimum {psi_min:10.5f} at ",end='')
    print(f"z = {z_min[0]:10.5f} and y = {y_min[0]:10.5f}.\n")
    
    extreme_loc = [(z_max,y_max),(z_min,y_min)]
    
# -- plot

    define_wall_shape( z*5.2, casecode=casecode, yshift=(h_ridge-h_md) )
    
    cbar = r'$\Psi$'
#    cbar_levels = np.linspace(-0.40,0.40,26)
#    cbar_ticks  = np.linspace(-0.4,0.4,5)

    for pure_stat in [True,False]:
        plot_slicex_stat( zz, yy, psi,
                        tag=tag,
                        filename=casecode+'_psi',
                        sonic=False,
                        wall=True,
                        col_map='RdBu_r',
                        cbar_label=cbar,
    #                      cbar_levels=cbar_levels,
    #                      cbar_ticks=cbar_ticks,
                        title=title,
                        extreme_loc=extreme_loc,
                        pure=pure_stat)

        plot_slicex_stat( zz, yy, w_favre,
                        tag=tag,
                        filename=casecode+'_w_favre',
                        sonic=False,
                        wall=True,
                        col_map='RdBu_r',
                        cbar_label=cbar,
    #                      cbar_levels=cbar_levels,
    #                      cbar_ticks=cbar_ticks,
                        title=title,
                        extreme_loc=extreme_loc,
                        pure=pure_stat)
        
        plot_slicex_stat( zz, yy, v_favre,
                        tag=tag,
                        filename=casecode+'_v_favre',
                        sonic=False,
                        wall=True,
                        col_map='RdBu_r',
                        cbar_label=cbar,
    #                      cbar_levels=cbar_levels,
    #                      cbar_ticks=cbar_ticks,
                        title=title,
                        extreme_loc=extreme_loc,
                        pure=pure_stat)


print(f"Finished doing slicing at x = {loc_delta:10.2f} delta.") 
print(f"=====================================\n")

# print out the time finishing the job
    
print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")