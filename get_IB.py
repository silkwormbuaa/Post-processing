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

from   vista.timer       import timer 

from   IB_force          import GetIBForce

from   ReadIn            import ReadZonegrpWall

#%% Get forces and pressure on IB
ForceFolderPath = "/media/wencanwu/Seagate Expansion Drive/temp/221221/forces_allblocks/"
ZonegrpFile = "/media/wencanwu/Seagate Expansion Drive/temp/221221/zonelist.dat"

os.chdir(ForceFolderPath)
os.chdir(os.pardir)
zonegrp = ReadZonegrpWall(ZonegrpFile)

with timer("Get forces on IB:"):
    GetIBForce(zonegrp,ForceFolderPath,"Cf_points_1221.dat")




#with timer("read one files") :
#    print(ReadAForce("f_bl_000794.dat"))
