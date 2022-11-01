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

import gc

from vista_tools         import *

from vista_pytecio       import *

from plt2pandas          import *

from timer               import timer


folder = "/home/wencanwu/my_simulation/temp/221014_lowRe/TP_stat"

outfolder = "/home/wencanwu/my_simulation/temp/221014_lowRe/stat_per_ave/"

listfile = "/home/wencanwu/my_simulation/temp/221014_lowRe/zonelist_sorted.dat"

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
                                   10.4,
                                   10.4,
                                   var_list)
        
        frame2tec3d(df,outfolder,str(i),zname=i)
        
        del df
        
        gc.collect()
        print("finish writing zone group %s"%i)
