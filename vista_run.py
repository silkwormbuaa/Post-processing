#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run.py
@Time    :   2022/10/15 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   executive script
'''

from vista_tools         import *
from vista_pytecio       import *
from timer               import timer


folderpath = "/home/wencanwu/my_simulation/temp/220825_lowRe/TP_stat"

listfile = "/home/wencanwu/my_simulation/temp/220825_lowRe/zonelist_sorted.dat"

with timer("get sorted zone list"):
    
#    get_zonegrp(folderpath)

    zonegrp = read_zonelist( listfile )
    
    regionrange = [-71.75, -1.2576, -60.6875, 22.8044]
    
    var_list = ['x','y','z','<u>','walldist']
       
    ave_block_ib( folderpath, 
                  zonegrp, 
                  "mean_result_test.dat", 
                  regionrange,
                  var_list)

