#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run.py
@Time    :   2022/10/15 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   get sorted zonelist and mean_result_ib_superficial average
'''

import sys

sys.path.append('..')

from utils.timer         import timer

from vista.pytecio       import *



folderpath = "/home/wencanwu/my_simulation/temp/220825_lowRe/TP_stat"

listfile = "/home/wencanwu/my_simulation/temp/220825_lowRe/zonelist_sorted.dat"

with timer("whole processing "):
    
    # output zonelist_sorted.dat
    
#    get_zonegrp(folderpath)

    zonegrps = read_zonelist( listfile )
    
    regionrange = [-57, -1.2576, -49.625, 22.8044]
    
    var_list = ['x',      'y',      'z',      '<u>',
                '<u`u`>', '<v`v`>', '<w`w`>', '<u`v`>',
                '<rho>' , 'walldist','<T>'              ]
       
    ave_block_ib( folderpath, 
                  zonegrps, 
                  "mean_result_ib_spf.dat", 
                  regionrange,
                  var_list,
                  spf_ave=True)

