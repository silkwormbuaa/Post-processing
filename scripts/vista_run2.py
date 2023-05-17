# -*- coding: utf-8 -*-
'''
@File    :   vista_run2.py
@Time    :   2022/10/19 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   script doing spanwise periodic averaging
'''

import os

import sys

source_dir = os.path.dirname( os.path.dirname( os.path.realpath(__file__) ))
sys.path.append( source_dir )

import gc

import tecplot             as     tp

from   vista.tools         import get_filelist

from   vista.timer         import timer

from   vista.pytecio       import read_zonelist

from   vista.pytecio       import spanwise_periodic_ave

from   plt2pandas          import frame2tec3d




folder = "/media/wencanwu/Seagate Expansion Drive/temp/221221/TP_stat"

outfolder = "/media/wencanwu/Seagate Expansion Drive/temp/221221/stat_per_ave/"

listfile = "/media/wencanwu/Seagate Expansion Drive/temp/221221/zonelist_sorted.dat"

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
