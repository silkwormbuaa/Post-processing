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

import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0] 
sys.path.append( source_dir )

from   vista.timer         import timer

from   vista.pytecio       import get_xline

from   vista.pytecio       import read_zonelist



#folderpath = "/home/wencanwu/my_simulation/temp/221125_lowRe/TP_stat"

#listfile = "/home/wencanwu/my_simulation/temp/221125_lowRe/zonelist_sorted.dat"

folderpath = "/media/wencanwu/Seagate Expansion Drive/temp/221221/TP_stat"

listfile = "/media/wencanwu/Seagate Expansion Drive/temp/221221/zonelist_sorted.dat"

with timer("whole processing "):
    
#    get_zonegrp(folderpath)

    zonegrps = read_zonelist( listfile )
    
    segment = [-0.260, -110.0, 115.0]
    
    var_list = ['x','y','z', 'walldist', '<p`p`>']
    
    get_xline( folderpath, zonegrps,
               segment, 
               "streamwise_line_1221.dat",
               var_list             )