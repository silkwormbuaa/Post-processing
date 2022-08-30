# -*- coding: utf-8 -*-
'''
@File    :   getCf.py
@Time    :   2022/08/30 16:50:49
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''



import tecplot as tp
import os
import numpy as np
from timer import timer

FilePath = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/"

os.chdir(FilePath)

dataset  = tp.data.load_tecplot("luis_1.dat")

x        = dataset.zone('ZONE 001').values('x').as_numpy_array()
y        = dataset.zone('ZONE 001').values('y').as_numpy_array()
mu       = dataset.zone('ZONE 001').values('<mu>').as_numpy_array()
u        = dataset.zone('ZONE 001').values('<u>').as_numpy_array()
#zones    = dataset.zone()
d_p      = 0.5*0.9886*507*507
Cf       = []
tau      = np.multiply(mu,np.divide(u,y))/d_p

with open("Cf.dat",'w') as f:
    for i in range(len(x)):
        f.write(str('%.7f'%x[i])+' '+str('%.7f'%tau[i])+'\n')
print(tau)

