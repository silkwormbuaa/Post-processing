# -*- coding: utf-8 -*-
'''
@File    :   run.py
@Time    :   2022/10/02 09:41:19
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Used as executable file.
'''
import os
import tecplot           as     tp
from   utils.timer       import timer
from   plt2pandas        import ReadPltBlocks_xy
from   plt2pandas        import frame2tec3d
from   ReadIn            import ReadZonegrp

Folder     = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/datafolder/'
SaveFolder = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/savefolder/'
pltfile    = '/home/wencanwu/my_simulation/temp/Low_Re_Luis/datafolder/datawithz.plt'

#with timer("read in file into tecplot dataset"):
#    dataset  = tp.data.load_tecplot(pltfile)

with timer("read in zonelist"):
    zonegrp = ReadZonegrp(Folder,'zonelist.dat')

with timer("read in dataset with tecplot"):
    dataset = tp.data.load_tecplot(pltfile)

with timer("read in plt into dataframe"):
    for i in range(len(zonegrp)):
        df = ReadPltBlocks_xy(dataset, zonegrp[i], SpanwiseAve=True)
        
        tecfilename = 'TP_stat_' + '{:06}'.format(i+1)
        frame2tec3d(df,SaveFolder,tecfilename,zname=(i+1))
#        print(df.index)
#        print(df['x'])
#        print(df['y'])
#        os.system("read -p 'Press Enter to continue...' var")
#with timer("remove overlapped points "):
#    df.drop_duplicates(subset=['x','y','z'],keep='first',inplace=True)
    
#with timer("reset index of df"):
#    df.reset_index(drop=True)
# use df.drop_duplicate(['x'],['y'])

#print(df['x'])
#print(df['y'])

#os.system("read -p 'Press Enter to continue...' var")


