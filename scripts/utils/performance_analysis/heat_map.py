#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   heat_map.py
@Time    :   2024/02/26 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import numpy             as     np
import pandas            as     pd
import seaborn           as     sb
import matplotlib.pyplot as     plt

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.tools       import get_filelist

def main():

    #display all rows
    pd.set_option('display.max_rows', None)
    
    out_dir = '/home/wencanwu/performance_comparison/roughwall_lowRe/1b1p/inca_out'
    
    files = get_filelist( out_dir, 'inca_0' )
    
    df_mean = get_mean( files )
    
    names = df_mean['Name'].tolist()
    
    # Inpsect output files
    # --------------------
    dfs = []
    for i, file in enumerate( files ):
        with open( file ) as f:
            
            content = [ line.rstrip('\n').split( '|' ) for line in ( f.readlines() ) if len( line ) > 100 and line[28] == '|' and line[1:5] != 'Name' ]
            
            df = pd.DataFrame( content, columns=['Name', 'Freq', 't_min', 't_avg', 't_max', 't_tot', 't_tot_pct'] )
            df = reset_order(df, names)
            df['Freq'] = df['Freq'].astype('int')
            for var in df.columns[2:]:
                df[var] = df[var].astype('float')
            dfs.append( df )

    # compose the data across different dataframes
        
    data = np.array([ np.array(df['t_tot_pct']) for df in dfs ])
    
    print(data.shape)
    print(len(names))
    
#    for i in range(len(names)):
#        data[:,i] = data[:,i]/mean_tavg[i]
    
    df = pd.DataFrame(data, columns=names, index=[ f'P{i}' for i in range(len(dfs))] )
    
    print(df)
    
    
    plt.figure(figsize=(32, 10))  
    sb.heatmap(df.T, cmap='coolwarm', fmt='.1f')  # annot=True 在热图上显示数值，cmap 设置颜色映射，fmt 格式化显示数字
    plt.title('Heatmap of DataFrame')  
    plt.subplots_adjust(bottom=0.3) 
    plt.show()  
    

# =============================================================================

def get_mean( files ):

    content = []
    for file in files:
       with open( file ) as f:
           content.extend( [ line.rstrip('\n').split( '|' ) for line in ( f.readlines() ) if len( line ) > 100 and line[28] == '|' and line[1:5] != 'Name' ] )

    # Reorder list alphabetically
    content.sort()
    
    df = pd.DataFrame( content, columns=['Name', 'Freq', 't_min', 't_avg', 't_max', 't_tot', 't_tot_pct'] )
   
    # convert to float
    df['Freq'] = df['Freq'].astype('int')
    for var in df.columns[2:]:
        df[var] = df[var].astype('float')
    
    df = df.groupby('Name').mean().reset_index()
    
    df.sort_values(by=['t_tot_pct'], ascending=False, inplace=True)
    
    return df

# =============================================================================

def reset_order(df:pd.DataFrame, name_order):
    
    df.set_index('Name', inplace=True)
    df_sorted = df.reindex(name_order, fill_value=0)
    df_sorted.reset_index(inplace=True)
    
    return df_sorted
    

# =============================================================================
# main program
# =============================================================================

if __name__ == "__main__": sys.exit( main() )