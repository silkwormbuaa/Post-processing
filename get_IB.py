# -*- coding: utf-8 -*-
'''
@File    :   get_IB.py
@Time    :   2022/09/25 23:01:04
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
from   timer             import timer 
from   IB_force          import *

#%% Get forces and pressure on IB
ForceFolderPath = "/home/wencanwu/my_simulation/temp/220825_lowRe/forces_allblocks/"
ZonegrpFile = "/home/wencanwu/my_simulation/temp/220825_lowRe/zonelist.dat"

os.chdir(ForceFolderPath)
os.chdir(os.pardir)
zonegrp = ReadZonegrpWall(ZonegrpFile)

with timer("Get forces on IB:"):
    GetIBForce(zonegrp,ForceFolderPath,"Cf_points.dat")




#with timer("read one files") :
#    print(ReadAForce("f_bl_000794.dat"))
