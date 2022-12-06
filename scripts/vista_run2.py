#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run2.py
@Time    :   2022/10/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   script doing spanwise periodic averaging
'''

import sys

sys.path.append('..')

import gc

from   utils.tools         import *

from   utils.timer         import timer

from   vista.pytecio       import *

from   plt2pandas          import *




folder = "/home/wencanwu/my_simulation/temp/221125_lowRe/TP_stat"

outfolder = "/home/wencanwu/my_simulation/temp/221125_lowRe/stat_per_ave/"

listfile = "/home/wencanwu/my_simulation/temp/221125_lowRe/zonelist_sorted.dat"

with timer("whole processing "):
    
#    get_zonegrp(folderpath)

    zonegrps = read_zonelist( listfile )
    
    var_list = ['x',      'y',      'z',      '<u>',
                '<v>',    '<w>',    '<T>',    '<rho>',
                '<u`u`>', '<v`v`>', '<w`w`>', '<u`v`>',
                '<p`p`>', '<p>',    'walldist'         ]


    
    filelist = sorted( get_filelist(folder) )
    
    print("Number of files in this folder:%s" % len(filelist))
    
    dataset = tp.data.load_tecplot_szl(filelist)

    os.chdir(outfolder)
    
    for i in range(len(zonegrps)):
        
#        dataset = tp.data.load_tecplot_szl(filelist)
        
        df = spanwise_periodic_ave(dataset,
                                   zonegrps[i],
                                   0.65,
                                   10.4,
                                   var_list)
        
        frame2tec3d(df,outfolder,str(i),zname=i)
        
        del df
        
        gc.collect()
        print("finish writing zone group %s"%i)
