#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   analyze_inca_performance.py
@Time    :   2024/01/22 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   script of reading outputfile of inspect_inca_out_files.py
             and do analysis
'''


import os
import sys
import pandas            as     pd

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

# ----------------------------------------------------------------------
# >>> main program
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/22  - created
#
# Desc
#
# ----------------------------------------------------------------------

def main():
    
    os.chdir('/home/wencanwu/Downloads/performance_compare/')
    
    file1 = '/home/wencanwu/Downloads/performance_compare/inca_out1/output1.out'
    file2 = '/home/wencanwu/Downloads/performance_compare/inca_out2/output2.out'
    
    df1 = read_inspect_output( file1 )
    df2 = read_inspect_output( file2 )
    
    # sort dataframe by name 
    
    df1.sort_values(by=['Name'], inplace=True)
    df2.sort_values(by=['Name'], inplace=True)
    
    pd.set_option('display.max_rows', None)
    
    # get data lists
    
    names = df1['Name'].tolist()
    
    freq1 = df1['Freq'].tolist()
    tmin1 = df1['t_min'].tolist()
    tmax1 = df1['t_max'].tolist()
    
    freq2 = df2['Freq'].tolist()
    tmin2 = df2['t_min'].tolist()
    tmax2 = df2['t_max'].tolist()
    
    r_tmax_21 = [tmax2[i]/tmax1[i] for i in range(len(df1))]
    
    data = {'Name': names, 'freq1':freq1, 'tmin1':tmin1, 'tmax1':tmax1, 
            'freq2':freq2, 'tmin2':tmin2, 'tmax2':tmax2, 'r_tmax_21':r_tmax_21}
    
    df = pd.DataFrame(data)
    
    # sort values and drop index
    
    df.sort_values(by=['r_tmax_21'], inplace=True, ascending=False, ignore_index=True)
    
    # display all rows
    
    print(df)
    
    df.to_string('inca_performance_compare.txt', index=False)

# ----------------------------------------------------------------------
# >>> read in inca_performance.out
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/22  - created
#
# Desc
#
# ----------------------------------------------------------------------

def read_inspect_output( filename ):
    
    with open(filename, 'r') as f:  lines = f.readlines()[5:]
        
    # read lines
    
    header = ['Name', 'Freq', 't_min', 't_max', 'time_ratio']
    
    df = pd.DataFrame(columns = header)
    
    for line in lines:
        
        data = line.strip().strip('\n').split()
        
        if data:  df.loc[len(df.index)] = data
                
    # convert data type
    
    df['Freq'] = df['Freq'].astype(int)
    df['t_min']     = df['t_min'].astype(float)
    df['t_max']     = df['t_max'].astype(float)
    df['time_ratio']= df['time_ratio'].astype(float)
    
    return df

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/22  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == '__main__': sys.exit( main() )