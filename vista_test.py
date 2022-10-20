#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_test.py
@Time    :   2022/10/20 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   For doing testing
'''
import os 

import tecplot           as      tp

import numpy             as      np 

import pandas            as      pd

from   vista_tools       import  *

from   vista_pytecio     import  *

from   timer             import  timer 

folder = "/home/wencanwu/my_simulation/temp/220927_lowRe/TP_stat"

outfolder = "/home/wencanwu/my_simulation/temp/220927_lowRe/stat_per_ave/"

listfile = "/home/wencanwu/my_simulation/temp/220927_lowRe/zonelist_sorted.dat"

with timer("whole processing "):
    
#    get_zonegrp(folderpath)

    zonegrps = read_zonelist( listfile )
    
    var_list = ['x',      'y',      'z',      '<u>',
                '<u`u`>', '<v`v`>', '<w`w`>', '<u`v`>',
                '<rho>' , '<p`p`>', 'walldist'         ]

    filelist = sorted( get_filelist(folder) )
    
    print("Number of files in this folder:%s" % len(filelist))
    
    dataset = tp.data.load_tecplot_szl(filelist)
    
    zonename1 =  zonegrps[0].zonelist[0][1]
    
    zonename2 =  zonegrps[0].zonelist[1][1]
    
    print(zonename1)
    print(zonename2)
    
    var = dataset.variable('z')
    var2 = dataset.variable('x')
    
    
    z = var.values(zonename1).as_numpy_array()
    x = var2.values(zonename1).as_numpy_array()
    
    z2 = var.values(zonename2).as_numpy_array()
    x2 = var2.values(zonename2).as_numpy_array()
    os.chdir(outfolder)
    
    np.savetxt("z.dat", z, delimiter=' ')    
    
    print(len(z))
    
    z = np.unique(z)
    z2 = np.unique(z2)
    
    z = np.column_stack((z,z2))
 
    x = np.unique(x)
    x2 =  np.unique(x2)
    x = np.column_stack((x,x2))
    
    np.savetxt("z_unique.dat",z ,delimiter=" ")
    
    np.savetxt("x_unique.dat",x ,delimiter=" ")

    print("length of z unique: %s"%len(z))
    
