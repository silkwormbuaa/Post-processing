#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run3.py
@Time    :   2022/10/24 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Script for getting a streamwise line data
'''


from vista_tools         import *

from vista_pytecio       import *

from timer               import timer


folderpath = "/home/wencanwu/my_simulation/temp/220927_lowRe/TP_stat"

listfile = "/home/wencanwu/my_simulation/temp/220927_lowRe/zonelist_sorted.dat"

with timer("whole processing "):
    
#    get_zonegrp(folderpath)

    zonegrps = read_zonelist( listfile )
    
    segment = [-0.468, -110.0, 115.0]
    
    var_list = ['x','y','z', 'walldist', '<p`p`>']
    
    get_xline( folderpath, zonegrps,
               -0.312,     segment, 
               "streamwise_line.dat",
               var_list             )