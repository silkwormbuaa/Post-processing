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

import sys

sys.path.append('..')

from utils.tools         import *

from utils.timer         import timer

from vista.pytecio       import *



folderpath = "/home/wencanwu/my_simulation/temp/221014_lowRe/TP_stat"

listfile = "/home/wencanwu/my_simulation/temp/221014_lowRe/zonelist_sorted.dat"

with timer("whole processing "):
    
#    get_zonegrp(folderpath)

    zonegrps = read_zonelist( listfile )
    
    segment = [-0.494, -110.0, 115.0]
    
    var_list = ['x','y','z', 'walldist', '<p`p`>']
    
    get_xline( folderpath, zonegrps,
               segment, 
               "streamwise_line_1014.dat",
               var_list             )