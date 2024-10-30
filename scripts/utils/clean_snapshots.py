#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   clean_snapshots.py
@Time    :   2024/10/30 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   clean snapshots with even index in snapshot_index.dat
'''

import os
import pandas            as    pd
import numpy             as    np

snapshots_dir = '/home/wencanwu/temp/241018/snapshots'

os.chdir(snapshots_dir)
index_file = 'snapshot_index.dat'

# read csv without header, but assgin the header manually
df = pd.read_csv(index_file, delimiter=r'\s+', header=None)

df.columns = ['index', 'time']

index = np.array(df['index'])

for i in index[1::2]:
    
    snap_folder = f'snapshot_{i:08d}'
    os.remove(snap_folder)
    print( f'{snap_folder} is removed.' )