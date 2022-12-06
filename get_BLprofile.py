# -*- coding: utf-8 -*-
"""
Original code credit to Weibo
"""

from   utils.timer       import timer

from   ReadIn            import GetZonegrp

from   ReadIn            import ReadZonegrp

from   BLprofile         import GetLine

from   vista.pytecio     import ave_block

from   vista.pytecio     import timeave_fbl

#%% Read plt data from INCA

FoldPath = "/home/wencanwu/my_simulation/temp/221125_lowRe/TP_stat"
#OutPath  = "/home/wencanwu/my_simulation/temp/220825_lowRe/DataPost/"
ForceFoldPath = "/home/wencanwu/my_simulation/temp/221125_lowRe/forces_3/"
with timer("all works"):
#    GetZonegrp(FoldPath)
    zonegrp = ReadZonegrp(FoldPath,'zonelist.dat')
#    ave_block(zonegrp,FoldPath,"mean_result_test.dat",3)
#    ReadBlock(zonegrp,FoldPath,"mean_result_test2.dat",3)
    timeave_fbl(ForceFoldPath)

#%% unfortunately, after loading one case, the dataset cannot
#   be released, so when trying to load another case, the 
#   loaded zones will go crazy.
"""
FoldPath = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/TP_stat"
OutPath  = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/DataPost/"
zonegrp2 = GetZonegrp(FoldPath)
with timer("get lines"):
    line_loc = [-58.25, -0.52, -58.25, 30.0]    
    GetLine(line_loc,zonegrp2,FoldPath,OutPath,1)
"""

'''
FoldPath = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/TP_stat"
OutPath  = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/DataPost/"

#zonegrp = GetZonegrp(FoldPath)
zonegrp = ReadZonegrp(FoldPath,'zonelist.dat')
with timer("get lines"):
    
    line_loc = [-30.0, 0.0, -30.0, 30.0]    
    GetLine(line_loc,zonegrp,FoldPath,OutPath,2)

    line_loc = [0.0, 0.0, 0.0, 30.0]    
    GetLine(line_loc,zonegrp,FoldPath,OutPath,2)

    line_loc = [50.4, 0.0, 50.4, 30.0]    
    GetLine(line_loc,zonegrp,FoldPath,OutPath,2)

    line_loc = [100.0, 0.0, 100.0, 30.0]    
    GetLine(line_loc,zonegrp,FoldPath,OutPath,2)

#    line_loc = [-71.7500, 0.0, -71.7500, 30.0]    
#    GetLine(line_loc,zonegrp,FoldPath,OutPath,2)
#
#    line_loc = [-68.0625, 0.0, -68.0625, 30.0]    
#    GetLine(line_loc,zonegrp,FoldPath,OutPath,2)
#
#    line_loc = [-64.3750, 0.0, -64.3750, 30.0]    
#    GetLine(line_loc,zonegrp,FoldPath,OutPath,2)
'''