#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   analyze_inca_performance.py
@Time    :   2024/01/22 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   script of reading inca_performance.out and analyze the performance
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
    
    file1 = '/home/wencanwu/Downloads/performance_compare/inca_performance1.out'
    file2 = '/home/wencanwu/Downloads/performance_compare/inca_performance2.out'
    
    df1 = read_inca_performance( file1 )
    df2 = read_inca_performance( file2 )
    
    # compare the 't_avg' of these two dataframes and get the iloc


#    index = [i for i in range(len(df1)) if (df1['t_avg'][i] < df2['t_avg'][i])]
    
    index = range(len(df1))
    
    names   = df1['Name'][index].tolist()
    t_avg1  = df1['t_avg'][index].tolist()
    t_freq1 = df1['Frequency'][index].tolist()
    t_avg2  = df2['t_avg'][index].tolist() 
    t_freq2 = df2['Frequency'][index].tolist()
    ratio   = [t_avg2[i]/t_avg1[i] for i in range(len(index))]
    time1   = [t_avg1[i]*t_freq1[i] for i in range(len(index))]
    time2   = [t_avg2[i]*t_freq2[i] for i in range(len(index))]
    
    data = {'Name': names, 't_avg1': t_avg1, 't_freq1': t_freq1, 't_avg2': t_avg2, 
            't_freq2': t_freq2, 'ratio': ratio, 'time': time1, 'time2': time2}
    
    df = pd.DataFrame(data)
    
    print(df)
    
    df.sort_values(by=['ratio'], inplace=True, ascending=False)
    
    # display all rows
    pd.set_option('display.max_rows', None)
    
    print(df)

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

def read_inca_performance( filename ):
    
    data_lines = []
    
    with open(filename, 'r') as f:
        
        lines = f.readlines()
        
        # if there are two blank lines, stop reading
        
        blank_line_count = 0
    
        # read lines
        
        for line in lines:
            
            line = line.strip().strip('\n')
            
            if line=='':
                blank_line_count += 1
                if blank_line_count == 2:
                    break
                
            else:
                
                # skip separator line
                if '-+-' not in line and 'Name' not in line:
                    data_lines.append( line )
                blank_line_count = 0
                
                
    # convert the list of lines to a string
    
    header = ['Name', 'Frequency', 't_min', 't_avg', 't_max', 't_tot', 't_tot_percent']
    
    df = pd.DataFrame(columns = header)
    
    for line in data_lines:
        line_content = [item.strip() for item in line.split('|')]
        df.loc[len(df.index)] = line_content
    
    # convert data type
    
    df['Frequency'] = df['Frequency'].astype(int)
    df['t_min']     = df['t_min'].astype(float)
    df['t_avg']     = df['t_avg'].astype(float)
    df['t_max']     = df['t_max'].astype(float)
    
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