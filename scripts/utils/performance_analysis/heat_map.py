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
from   matplotlib.colors import LogNorm, Normalize

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.tools       import get_filelist

# ----------------------------------------------------------------------
# >>> evaluate_local_performance_fluctuation                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def evaluate_local_performance_fluctuation( out_dir ):
    
    """
    output the heatmap of tmax/tmin for each processor each subroutine.
    
    fluctuation may result from the following reasons:
    1. other subroutines' performance fluctuation may affect the current subroutine's performance,
       e.g. execute immediately or wait for other processors.
    2. some tasks within this subroutine may/may not be executed. Like post-processing.
    
    """
    
    files = get_filelist( out_dir, 'inca_0' )
    
    df_mean = get_mean( files )          # df are sorted based on t_tot_pct
    df_mean.reset_index(drop=True, inplace=True)
    
    names = df_mean['Name'].tolist()
    
    print(names)
    pd.set_option('display.max_rows', None)
    print( df_mean )

# ---- read in inca.out files from all processors

    dfs = read_incaouts2df( files, names )    

# ----- compose t_max/t_min data frame
        
    data = np.array([ np.array(df['t_max'])/(np.array(df['t_min'])+0.00000001) for df in dfs ])    
    df_t_tot_pct = pd.DataFrame(data, columns=names, index=[ f'P{i:04d}' for i in range(len(dfs))] )
    df_t_tot_pct = df_t_tot_pct.iloc[:,1:33]
    
    plt.figure(figsize=(64, 10))  
    sb.heatmap(df_t_tot_pct.T, cmap='coolwarm', fmt='.1f', norm=LogNorm(1,100000))
    plt.title('Heatmap of local t_max/t_min')  
    plt.subplots_adjust(bottom=0.3) 
    plt.savefig('local_performance_fluctuation.png')
    plt.show()
    plt.close()
    
    
# ----------------------------------------------------------------------
# >>> compare_2incaout                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def compare_2incaout( out_dir, out_dir2):
    
    """
    show the heatmap of t_avg_2/t_avg_1 for each processor each subroutine.
    """

    files = get_filelist( out_dir, 'inca_0' )
    files2 = get_filelist( out_dir2, 'inca_0' )
    
    df_mean = get_mean( files2 )          # df are sorted based on t_tot_pct
    df_mean.reset_index(drop=True, inplace=True)
    
    names = df_mean['Name'].tolist()
    
    print(names)
    pd.set_option('display.max_rows', None)
    print( df_mean )
    
# ---- read in inca.out files from all processors

    dfs = read_incaouts2df( files, names )
    dfs2 = read_incaouts2df( files2, names )
    
# ----- compare each processor's performance

    data = np.array([ np.array(df2['t_avg']) / (np.array(df['t_avg'])+0.0000000001) for df,df2 in zip(dfs,dfs2) ])
    
    df_2t_avg = pd.DataFrame(data, columns=names, index=[ f'P{i:04d}' for i in range(len(dfs))] )
    df_2t_avg = df_2t_avg.iloc[:,0:33]
    
    plt.figure(figsize=(64, 10))  
    sb.heatmap(df_2t_avg.T, cmap='Greys', fmt='.1f',norm=LogNorm(1,100))
    plt.title('Heatmap of 2t_avg')  
    plt.subplots_adjust(bottom=0.3) 
    plt.savefig('2t_avg.png')
    plt.show()
    plt.close()
    
# ----------------------------------------------------------------------
# >>> show_1incaout                                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def show_1incaout( out_dir ):

    """
    show the total time percentage distribution from all processors.
    show the average time normalized by the average value of corresponding subroutine. 
    """
    os.chdir( out_dir )
    
    files = get_filelist( out_dir, 'inca_0' )
    
    df_mean = get_mean( files )          # df are sorted based on t_tot_pct
    df_mean.reset_index(drop=True, inplace=True)
    
    names = df_mean['Name'].tolist()
    
    print(names)
    pd.set_option('display.max_rows', None)
    print( df_mean )
    
# ---- read in inca.out files from all processors

    dfs = read_incaouts2df( files, names )    

# ----- compose t_tot_pct data frame
        
    data = np.array([ np.array(df['t_tot_pct']) for df in dfs ])    
    df_t_tot_pct = pd.DataFrame(data, columns=names, index=[ f'P{i:04d}' for i in range(len(dfs))] )
    df_t_tot_pct = df_t_tot_pct.iloc[:,3:33]
    
    plt.figure(figsize=(64, 10))  
    sb.heatmap(df_t_tot_pct.T, cmap='coolwarm', fmt='.1f', vmin=0, vmax=50)
    plt.title('Heatmap of t_tot_pct')  
    plt.subplots_adjust(bottom=0.3) 
    plt.savefig('t_tot_pct.png')
    plt.show()
    plt.close()
    
# ----- compose t_avg data frame

    data = np.array([ np.array(df['t_avg']) for df in dfs ])
    
    # normalize t_avg by the average of t_avg
    
    for i in range(len(names)):
        data[:,i] = data[:,i]/df_mean['t_avg'][i]
    
    df_t_avg = pd.DataFrame(data, columns=names, index=[ f'P{i:04d}' for i in range(len(dfs))] )
    df_t_avg = df_t_avg.iloc[:,0:33]
    
    plt.figure(figsize=(64, 10))
    sb.heatmap(df_t_avg.T, cmap='coolwarm', fmt='.1f',vmin=0, vmax=5)
    plt.title('Heatmap of t_avg')
    plt.subplots_adjust(bottom=0.3)
    plt.savefig('t_avg.png')
    plt.show()
    plt.close()


# ----------------------------------------------------------------------
# >>> get_mean                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def get_mean( files ):
    
    """
    get the processor-wise average of all subroutines and time type.
    """

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

# ----------------------------------------------------------------------
# >>> reset_order                                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def reset_order(df:pd.DataFrame, name_order):
    
    """
    reset the order of the data frame based on the given name_order.
    """
    
    df.set_index('Name', inplace=True)
    df_sorted = df.reindex(name_order, fill_value=0)
    df_sorted.reset_index(inplace=True)
    
    return df_sorted


# ----------------------------------------------------------------------
# >>> read_incaouts2df                                              (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/13  - created
#
# Desc
#
# ----------------------------------------------------------------------

def read_incaouts2df( files, names ): 
    
    """
    read in all given inca.out files and reset data order based on the given names.
    
    return sorted data frames list.
    """

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
            
    return dfs


# =============================================================================
# main program
# =============================================================================

if __name__ == "__main__": 
    
    outpath = '/home/wencanwu/performance_comparison/midRe/'
    
    os.chdir( outpath )
    
    out_dir = '/home/wencanwu/performance_comparison/midRe/genoa/inca_out'
    out_dir2 = '/home/wencanwu/performance_comparison/midRe/rome/inca_out'
    
    # compare_2incaout( out_dir, out_dir2 )
    
    # show_1incaout( out_dir )
    # show_1incaout( out_dir2 )

    evaluate_local_performance_fluctuation( out_dir )
    # evaluate_local_performance_fluctuation( out_dir2 )
