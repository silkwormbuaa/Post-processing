#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   get_wd.py
@Time    :   2023/09/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


import os
import sys
import pickle
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt
import matplotlib.colors as     colors

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.statistic   import StatisticData

from   vista.snapshot    import Snapshot

from   vista.grid        import GridData

from   vista.timer       import timer

from   vista.plane_analy import save_isolines

from   vista.lib.form    import phy
from   vista.lib.form    import mth

snapshotfile = '/home/wencanwu/my_simulation/temp/220927_lowRe/wall_dist/snapshot_01329504/snapshot.bin'
resultspath  = '/home/wencanwu/my_simulation/temp/220927_lowRe/results'

stat_file = resultspath + '/statistics.bin'
grid_file = resultspath + '/inca_grid.bin'
ccfile = resultspath + '/cutcells_setup.dat'

bbox = [ -60.0, 100.0, -1.3, 0.01, -11.0, 11.0]

# read in statistics.bin and also grid and cut cell info.
with timer("read in grid"):
    G = GridData( grid_file )
    G.read_grid()
    
    block_list = G.select_blockgrids( bbox, mode='within' )

with timer("read statistics"):
    S = StatisticData( stat_file )

    vars = ['u','mu']

    with open( stat_file, 'br' ) as f:
        
        S.read_stat_header( f )
        S.read_stat_body(f, block_list, vars)
    
    S.grid3d = G
    
# read in wall distance

wd_snap = Snapshot( snapshotfile )

wd_snap.verbose = False

wd_snap.read_snapshot(block_list)

for num in block_list:
    
    S.bl[num-1].df['wd'] = wd_snap.snap_data[num-1][5]['wd']
    
print(S.bl[num-1].df)

with timer("read cutcell info"):
    
    cc_df = pd.read_csv( ccfile, delimiter=r'\s+')

    cc_df.drop( columns=['fax0','fax1','faz0','faz1','processor'], inplace=True )

with timer("compute friction projection"):
    
    S.friction_projection( block_list, G, cc_df )
    


xx = np.array( S.df_fric['x'] )
zz = np.array( S.df_fric['z'] )
fric = np.array( S.df_fric['fric'] )

npx = len( np.unique(xx) )
npz = len( np.unique(zz) )

xx = xx.reshape(npz,npx)
zz = zz.reshape(npz,npx)
fric = fric.reshape(npz,npx)

save_isolines(xx,zz,fric,0.0,"separation_lines.pkl")

fig, ax = plt.subplots()

cs = ax.contourf(xx,zz,fric,levels=51,cmap='coolwarm',norm=colors.CenteredNorm())

cbar = plt.colorbar(cs)

with open('separation_lines.pkl','rb') as f:
    lines = pickle.load( f )
    
    for line in lines:
        x_sep = line[:,0]
        z_sep = line[:,1]
        ax.plot(x_sep,z_sep,'black',linewidth=0.8)

plt.show()


    
#    for num in block_list:
#        
#        temp_df = cc_df[ cc_df['block_number'] == num ]
#        
#        print(temp_df)
#        
#        G.g[num-1].assign_vol_fra( temp_df, geo_case )
#
#with timer("get a df with same i,k"):
#    
##    for i in range(1,38):
##        for k in range(1,38):
##    
##            df = temp_df[ (temp_df['i'] == i) & (temp_df['k'] == k) ]
#    group = temp_df.groupby(['i','k'])
#
#    print(group.ngroups)
#    
#    for i in range(4,36):
#        for k in range(4,36):
#            try:
#                df = group.get_group((i,k))
#            except KeyError:
#                print(f"{i},{k} no value")

            