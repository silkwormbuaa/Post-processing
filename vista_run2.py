#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_run2.py
@Time    :   2022/10/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''


from vista_tools         import *

from vista_pytecio       import *

from plt2pandas          import *

from timer               import timer


folder = "/home/wencanwu/my_simulation/temp/220927_lowRe/TP_stat"

outfolder = "/home/wencanwu/my_simulation/temp/220927_lowRe/stat_per_ave/"

listfile = "/home/wencanwu/my_simulation/temp/220927_lowRe/zonelist_sorted.dat"

with timer("whole processing "):
    
#    get_zonegrp(folderpath)

    zonegrps = read_zonelist( listfile )
    
#    var_list = ['x',      'y',      'z',      '<u>',
#                '<u`u`>', '<v`v`>', '<w`w`>', '<u`v`>',
#                '<rho>' , '<p`p`>', 'walldist'         ]

    var_list = ['x','y','z','<u>']
    
    filelist = sorted( get_filelist(folder) )
    
    print("Number of files in this folder:%s" % len(filelist))
    
    dataset = tp.data.load_tecplot_szl(filelist)

#    for i in range(len(zonegrps)):
    os.chdir(outfolder)
    for i in range(1):    
        df = spanwise_periodic_ave(dataset,
                                   zonegrps[i],
                                   1.3,
                                   5.2,
                                   var_list)
        frame2tec3d(df,outfolder,str(i))
        print("finish writing zone group %s"%i)
