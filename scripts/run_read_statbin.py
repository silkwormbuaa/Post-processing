#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   read_statbin.py
@Time    :   2023/02/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   read statistic binary data.
           ! Set the headers for cutcell_setup first
'''

import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

import pandas            as     pd

import numpy             as     np

from   vista.statistic   import StatisticData

from   vista.grid        import GridData

from   vista.timer       import timer

#datapath = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results/'

datapath = '/home/wencanwu/my_simulation/temp/220927_lowRe/results/'

datafile = datapath + 'statistics.bin'
gridfile = datapath + 'inca_grid.bin'

outpath  = datapath

ccfile   = datapath + 'cutcells_setup.dat'

outfile  = 'mean_profile_test.dat'

# - select which wavy wall case
#
#   1 : 1014 case, D/delta = 2
#   2 : 0926 case, D/delta = 1
#   3 : 0825 case, D/delta = 0.5
#   4 : 0927 case, D/delta = 0.25
#   5 : 1221 case, D/delta = 0.125

geo_case = 4

G = GridData( gridfile )

#os.chdir( outpath )

G.verbose = True

## read in whole grid info

# 1. read_grid() : read_grid_header() + read_grid_body()
#   * grid headers: containing what will be read
#   * grid body: every block's grids information
# 2. sorted_group: sort and group block index basd on x,y,z

G.read_grid()

G.get_sorted_groups()


# given a rectangular region, get a list of blocks within the region

rect1 = [-57, -1.2576, -49.625, 1.73695342]

G.select_blockgrids( rect1 )

block_list = np.array( G.blockgrids_sel ).ravel()

with timer("read in cut cell info "):

    cc_df = pd.read_csv( ccfile, delimiter = r'\s+' )

    cc_df.drop( columns=['nx',  'ny',  'nz'
                        ,'fax0','fay0','faz0'
                        ,'fax1','fay1','faz1']
                        , inplace=True )
    
with timer("assign volume fractions "):
    
    for num in block_list:

        # dataframe slice for a certain block
        temp_df = cc_df[cc_df['block_number'] == num ]
        
        # block number starts from 1, but python list index
        # starts from 0
        
        G.g[num-1].assign_vol_fra( temp_df, geo_case )

S = StatisticData( datafile )

with timer("read block statistics data "):
    with open( datafile, 'br' ) as f:   
            
        S.read_stat_header( f )

        # only blocks in the list will be filled data chunk
        S.read_stat_body( f, block_list )
        

"""
with timer("calculating profile "):
    
    npy = G.g[block_list[0]-1].npy
    
    u_ls   = list()
    rho_ls = list()
    T_ls   = list()
    uu_ls  = list()
    uv_ls  = list()
    vv_ls  = list()
    ww_ls  = list()
    
    # ignore warning when numpy does 0/0 operation 
    np.seterr(invalid='ignore')
    
    for j in range(3,npy+3):
        
        vol = list()
        u   = list()
        rho = list()
        T   = list()
        uu  = list()
        uv  = list()
        vv  = list()
        ww  = list()
        
        # now the code can only do averaging in the 
        #     bottom block
        
        
        for num in block_list:
            
            vol.append( G.g[num-1].vol_fra[:,j,:] )
            u.append( S.bl[num-1].mean[:,j,:,0])
            rho.append( S.bl[num-1].mean[:,j,:,3])
            T.append( S.bl[num-1].mean[:,j,:,6])
            uu.append( S.bl[num-1].uu[:,j,:])
            uv.append( S.bl[num-1].uv[:,j,:])
            vv.append( S.bl[num-1].vv[:,j,:])
            ww.append( S.bl[num-1].ww[:,j,:])
            
        vol = np.array(vol)
        u   = np.array(u)
        rho = np.array(rho)
        T   = np.array(T)
        uu  = np.array(uu)
        uv  = np.array(uv)
        vv  = np.array(vv)
        ww  = np.array(ww)
        
        up = np.multiply(vol,u)        
        u_mean = np.divide(np.sum(up),np.sum(vol))
        u_mean = np.nan_to_num(u_mean)
        u_ls.append( u_mean )

        up = np.multiply(vol,rho)        
        rho_mean = np.divide(np.sum(up),np.sum(vol))
        rho_mean = np.nan_to_num(rho_mean)
        rho_ls.append( rho_mean )

        up = np.multiply(vol,T)        
        T_mean = np.divide(np.sum(up),np.sum(vol))
        T_mean = np.nan_to_num(T_mean)
        T_ls.append( T_mean )
        
        up = np.multiply(vol,uu)        
        uu_mean = np.divide(np.sum(up),np.sum(vol))
        uu_mean = np.nan_to_num(uu_mean)
        uu_ls.append( uu_mean )

        up = np.multiply(vol,uv)        
        uv_mean = np.divide(np.sum(up),np.sum(vol))
        uv_mean = np.nan_to_num(uv_mean)
        uv_ls.append( uv_mean )

        up = np.multiply(vol,vv)        
        vv_mean = np.divide(np.sum(up),np.sum(vol))
        vv_mean = np.nan_to_num(vv_mean)
        vv_ls.append( vv_mean )
        
        up = np.multiply(vol,ww)        
        ww_mean = np.divide(np.sum(up),np.sum(vol))
        ww_mean = np.nan_to_num(ww_mean)
        ww_ls.append( ww_mean )
#        u_profile = np.nan_to_num(u_profile)

    y_ls = G.g[block_list[0]-1].gy[3:npy+3]
    
os.chdir( outpath )
    
with open(outfile, 'w') as f:
    f.write('{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}{:<17}'\
            .format('y', '<u>', '<rho>', '<T>', \
            '<u`u`>','<u`v`>','<v`v`>' ,'<w`w`>') + '\n')
    
    for i in range(len(y_ls)):
        f.write(str('{:<17.8e}'.format(y_ls[i]))  )
        f.write(str('{:<17.8e}'.format(u_ls[i]))  )
        f.write(str('{:<17.8e}'.format(rho_ls[i])))
        f.write(str('{:<17.8e}'.format(T_ls[i]))  )
        f.write(str('{:<17.8e}'.format(uu_ls[i])) )
        f.write(str('{:<17.8e}'.format(uv_ls[i])) )
        f.write(str('{:<17.8e}'.format(vv_ls[i])) )
        f.write(str('{:<17.8e}'.format(ww_ls[i])) + '\n' )
        
"""
        