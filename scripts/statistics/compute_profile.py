#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   read_statbin.py
@Time    :   2023/02/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Extract profile from statistic binary data.
           ! Set the headers for cutcell_setup first
'''

import os
import sys
import time

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

import pandas            as     pd
import numpy             as     np

from   vista.log         import Logger
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.params      import Params
from   vista.colors      import colors          as    clr
from   vista.snapshot    import Snapshot
from   vista.statistic   import StatisticData
from   vista.directories import Directories
from   vista.directories import create_folder
from   vista.tools       import get_filelist
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================
# option zone
# =============================================================================

locs = [-110, -79.6, -53.6, -27.6, 76.4, 102.4] 

## locs in dimensionless [-31.33, -25, -20,-15, 5, 10 ] 

# =============================================================================

dirs = Directories( os.getcwd() )

outpath  = dirs.pp_profile_array
datafile = dirs.statistics
gridfile = dirs.grid
ccfile   = dirs.cc_setup
parametersfile = dirs.case_para_file

# read parameters to check if it is a rough wall case
params     = Params( parametersfile )
roughwall  = params.roughwall

if roughwall:
    snapshotfile = get_filelist(dirs.wall_dist,key='snapshot.bin')[0]

# - enter outpath

create_folder( outpath )
os.chdir(outpath)

# =============================================================================

for i, loc in enumerate(locs):
    
    print(f"Processing extracting profile at location: {loc}.\n")
    
    bbox    = [ loc-2.5, loc+2.5, -1.2576, 86.0,  -11.0, 11.0]  # bounding box
    wbox    = [ loc-2.5, loc+2.5, -1.2576, 0.0,   -11.0, 11.0]  # bounding box for rough wall
    wbox_sw = [ loc-2.5, loc+2.5, 0.0,     0.01,  -11.0, 11.0]

    outfile  = f'profile_mean_{i:02d}.dat'
    wall_stat_file = f'wall_statistics_{i:02d}.dat'

    # - read in grid file

    with timer("read in grid"):
        
        G = GridData( gridfile )
        G.verbose = False

        G.read_grid()
        G.cell_volume()
        
        # given a rectangular region, get a list of blocks in the region
        block_list = G.select_blockgrids( bbox )
        
        if roughwall:
            block_list_w = G.select_blockgrids( wbox )
        else: block_list_w = G.select_blockgrids( wbox_sw )


    # - read in statistics data

    with timer("read block statistics data "):

        S = StatisticData( datafile )
        
        with open( datafile, 'br' ) as f:   
                
            S.read_stat_header( f )
            vars = ['u','v','w','p','pp','rho','mu','T','uu','vv','ww','uv']
            
            # only blocks in the list will be filled data chunk
            S.read_stat_body( f, block_list, vars )


    # - assign vol_fra to grid, then match G to S

    with timer("Assign vol_fra to G"):
        
        if roughwall:

            # - read in wall distance data
            with timer("read wall distance field"):
                
                wd_snap = Snapshot( snapshotfile )
                wd_snap.read_snapshot( block_list, var_read=['wd'] )

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

        # match grid and pass vol_fra from G to S
        S.match_grid( block_list, G )


    # - compute wall friction

    with timer("compute wall friction"):
        
        if roughwall:

            # match wall distance to S
            for num in block_list:
                S.bl[num-1].df['wd'] = wd_snap.snap_data[num-1].df['wd']
                
            #print(S.bl[num-1].df)
            
            S.friction_projection( block_list_w, G, cc_df )
            S.wall_vars_projection( block_list_w, G, cc_df )

            S.df_fric = S.df_fric[ (S.df_fric['x']>wbox[0]) & (S.df_fric['x']<wbox[1]) ]
            S.df_wall = S.df_wall[ (S.df_wall['x']>wbox[0]) & (S.df_wall['x']<wbox[1]) ]
            
            tau_ave = np.array(S.df_fric['fric']).mean()
            mu_ave  = np.array(S.df_wall['mu']).mean()
            rho_ave = np.array(S.df_wall['rho']).mean()

            u_tau = np.sqrt( abs(tau_ave/rho_ave) )
            lv    = mu_ave/rho_ave/u_tau    
                
            with open(wall_stat_file,'w') as f:
                
                f.write(f"tau_ave  {tau_ave:15.5f}\n")
                f.write(f"rho_ave  {rho_ave:15.5f}\n")
                f.write(f"mu_ave   {mu_ave:15.5f} \n")
                f.write(f"u_tau    {u_tau:15.5f}  \n")
                f.write(f"lv       {lv:15.5f}     \n")
            
            # remove useless file generated by projection function
            os.remove("friction_projection.pkl")
            os.remove("wall_vars_projection.pkl")
            
        else:
            
            S.extract_wall_vars_sw( block_list_w, G )
            
            S.df_wall = S.df_wall[ (S.df_wall['x']>wbox[0]) & (S.df_wall['x']<wbox[1]) ]
            
            tau_ave = np.array(S.df_wall['fric']).mean()
            mu_ave  = np.array(S.df_wall['mu']).mean()
            rho_ave = np.array(S.df_wall['rho']).mean()

            u_tau = np.sqrt( abs(tau_ave/rho_ave) )
            lv    = mu_ave/rho_ave/u_tau    
                
            with open(wall_stat_file,'w') as f:
                
                f.write(f"tau_ave  {tau_ave:15.5f}\n")
                f.write(f"rho_ave  {rho_ave:15.5f}\n")
                f.write(f"mu_ave   {mu_ave:15.5f} \n")
                f.write(f"u_tau    {u_tau:15.5f}  \n")
                f.write(f"lv       {lv:15.5f}     \n")
            
            # remove useless file generated by projection function
            os.remove("wall_vars_projection.pkl")
            

    # - with vol_fra assigned, compute profile

    with timer("compute profile"):
        
        S.drop_ghost( block_list )
        S.compute_profile( block_list, bbox, vars, 
                           outfile=outfile, roughwall=roughwall )

    del S
    
    # - remind to modify the profile
    
    print(f"Saved {outfile} and {wall_stat_file}.")
    
    if not roughwall:
        print(clr.bg.red,
            "Please add the values at the wall(y=0) manually!",
            clr.reset)
    else:
        print(clr.bg.red,
            f"Please modify the wall coordinate(y=-{params.H}) manually!",
            clr.reset,'\n')

    print(f"Finished extracting profile at location {loc}.\n")
# print out the time finishing the job

print(f"Finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
sys.stdout.flush()    