#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   integrate_slice.py
@Time    :   2024/03/08 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Intergrate variable over outlet slice.
'''


import os
import sys
import pickle
import numpy             as     np

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.statistic   import StatisticData
from   vista.grid        import GridData
from   vista.directories import Directories
from   vista.timer       import timer
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

locs = [-115.0, 119.9 ]
slic_type = 'X'         # only 'X'
range = [-10.4, 0.0, 10.4, 10]

opt = 2 

# =============================================================================

dirs = Directories( os.getcwd() )

datafile = dirs.statistics
gridfile = dirs.grid
outpath  = dirs.pp_integral

# - read in grid info

G = GridData( gridfile)
G.read_grid()

# - enter output directory
if not os.path.exists( outpath ): os.mkdir( outpath )

os.chdir( outpath )

# - do slicing and integrate

for i, loc in enumerate( locs ):
    
    df_slice_file = f"df_slice_{slic_type}_{i:02d}.pkl"
    output_file   = f"intergral_{slic_type}_{i:02d}.dat"
    
    if not os.path.exists( df_slice_file ):
        
        print(f"Start doing slicing at {slic_type} = {loc:10.2f}.\n")
        
        block_list, indx_slic = G.select_sliced_blockgrids( slic_type, loc )
        
        print(f"Selected {len(block_list)} blocks.\n")
        
        # - read statistics data file
        
        S = StatisticData( datafile )
        
        with timer("Read statistics data file and match grid"):
            
            with open( datafile, 'rb' ) as f:
                
                S.read_stat_header( f )
                vars = ['u','v','w','rho','p','T','urho','uu','vv','ww']
                S.read_stat_body( f, block_list, vars )
                S.match_grid( block_list, G )
        
        with timer("Get slice dataframe"):
            
            df_slice =  S.get_slice_df( block_list, G, indx_slic, slic_type )
            
            with open( df_slice_file, 'wb' ) as f:
                pickle.dump( df_slice, f )
    
    # - read in slice dataframe
    else:
        print(f"{df_slice_file} already exists, read in directly...\n")
        df_slice = pickle.load( open( df_slice_file, 'rb' ) )

# ----- integrate

    with timer("Integrate"):
        
        print(df_slice)
        
        df_slice.drop(  df_slice[ (df_slice['z']<range[0]) |
                                  (df_slice['y']<range[1]) |
                                  (df_slice['z']>range[2]) |
                                  (df_slice['y']>range[3]) ].index, 
                        inplace=True )
        
        df_slice['area'] = np.array(df_slice['hy']) * np.array(df_slice['hz'])
        print( df_slice )
        
        R  = 287.0508571
        gamma = 1.4
        Cp = R*gamma/(gamma-1)
        
        if opt == 1:
            ke = 0.5*( np.array(df_slice['u'])**2 + 
                       np.array(df_slice['v'])**2 + 
                       np.array(df_slice['w'])**2)
            df_slice['Tt'] = np.array(df_slice['T']) + ke/Cp
            df_slice['pt'] = np.array(df_slice['p']) * \
                             (np.array(df_slice['Tt'])/np.array(df_slice['T']))**\
                             (gamma/(gamma-1))        
            df_slice['m'] = np.array(df_slice['area'])*np.array(df_slice['rho'])*np.array(df_slice['u'])
        
        elif opt == 2:
            ke = 0.5*(np.array(df_slice['uu']) + np.array(df_slice['vv']) + np.array(df_slice['ww']))
            df_slice['Tt'] = np.array(df_slice['T']) + ke/Cp
            df_slice['pt'] = np.array(df_slice['p']) * \
                             (np.array(df_slice['Tt'])/np.array(df_slice['T']))**\
                             (gamma/(gamma-1))
            df_slice['m'] = np.array(df_slice['area'])*np.array(df_slice['urho'])

        m_total = np.sum(df_slice['m'])
        pt_avg = np.sum(df_slice['pt'] * df_slice['m']) / m_total
        Tt_avg = np.sum(df_slice['Tt'] * df_slice['m']) / m_total
        
        print(f"m_total = {m_total:15.2f}")
        print(f"Tt_avg  = {Tt_avg:15.2f}")
        print(f"pt_avg  = {pt_avg:15.2f}")
        
        with open( output_file, "w" ) as f:
            f.write(f"m_total = {m_total:15.2f}\n")
            f.write(f"Tt_avg  = {Tt_avg:15.2f}\n")
            f.write(f"pt_avg  = {pt_avg:15.2f}\n")
