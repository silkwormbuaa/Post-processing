#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   txt.py
@Time    :   2022/11/29 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Methods used for manipulating texts.
'''


# ----------------------------------------------------------------------
# >>> Replace text                                              ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/29  - created
#
# Desc
# 
# - replace the text 
#
# ----------------------------------------------------------------------

def txt_replace( file, in_str, out_str ):

    with open( file, 'r' ) as f:
        
        lines = f.readlines()
        
        content = ''
        
        for line in lines:
            
            line = line.replace( in_str, out_str )
            
            content = content + line
        
        print( content )

    with open( file, 'w' ) as f:
        
        f.write( content )
        
#%% executable test scripts

if __name__=='__main__':
    
    pass