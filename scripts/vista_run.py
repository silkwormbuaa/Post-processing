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

import sys

sys.path.append('..')

from utils.tools         import *

from utils.timer         import timer

from vista.pytecio       import *



folderpath = "/home/wencanwu/my_simulation/temp/221014_lowRe/TP_stat"

listfile = "/home/wencanwu/my_simulation/temp/221014_lowRe/zonelist_sorted.dat"

with timer("whole processing "):
    
    get_zonegrp(folderpath)

    zonegrps = read_zonelist( listfile )
    
    regionrange = [-71.75, -1.2576, -60.6875, 22.8044]
    
    var_list = ['x',      'y',      'z',      '<u>',
                '<u`u`>', '<v`v`>', '<w`w`>', '<u`v`>',
                '<rho>' , 'walldist'                   ]
       
    ave_block_ib( folderpath, 
                  zonegrps, 
                  "mean_result_ib_spf.dat", 
                  regionrange,
                  var_list,
                  spf_ave=True)

