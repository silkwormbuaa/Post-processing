# -*- coding: utf-8 -*-
'''
@File    :   getCf.py
@Time    :   2022/08/30 16:50:49
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os

import tecplot           as     tp

import numpy             as     np

from   utils.timer       import timer

# from tecplot line to get Cf and Cp

FilePath = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/"

os.chdir(FilePath)

dataset  = tp.data.load_tecplot("luis_3.dat")

x        = dataset.zone('ZONE 001').values('x').as_numpy_array()
y        = dataset.zone('ZONE 001').values('y').as_numpy_array()
mu       = dataset.zone('ZONE 001').values('<mu>').as_numpy_array()
u        = dataset.zone('ZONE 001').values('<u>').as_numpy_array()
p        = dataset.zone('ZONE 001').values('<p>').as_numpy_array()
#zones   = dataset.zone()
d_p      = 0.5*0.9886*507*507
p_i      = 45447.289
Cf       = np.multiply(mu,np.divide(u,y))/d_p
Pf       = np.divide(p,p_i)
x        = (x-50.4)/7.67277

with open("Cf_flat_all.dat",'w') as f:
    for i in range(len(x)):
        f.write(str('%.7f'%x[i])+' '+str('%.7f'%Cf[i])+' '\
              + str('%.7f'%Pf[i])+'\n')
#print(tau)
'''
FilePath = "/home/wencanwu/my_simulation/temp/220825_lowRe/linedata/"
#---before run this code, change zonename in data file first.
os.chdir(FilePath)

dataset   = tp.data.load_tecplot("crest_w.dat")
dataset   = tp.data.load_tecplot("crest_2.dat")
x1        = dataset.zone('ZONE001').values('x').as_numpy_array()
y1        = dataset.zone('ZONE001').values('y').as_numpy_array()
mu1       = dataset.zone('ZONE001').values('<mu>').as_numpy_array()
u1        = dataset.zone('ZONE001').values('<u>').as_numpy_array()
p1        = dataset.zone('ZONE001').values('<p>').as_numpy_array()

x2        = dataset.zone('ZONE002').values('x').as_numpy_array()
y2        = dataset.zone('ZONE002').values('y').as_numpy_array()
mu2       = dataset.zone('ZONE002').values('<mu>').as_numpy_array()
u2        = dataset.zone('ZONE002').values('<u>').as_numpy_array()
p2        = dataset.zone('ZONE002').values('<p>').as_numpy_array()
#zones    = dataset.zone()
d_p       = 0.5*0.9886*507*507
p_i       = 45447.289
du        = np.subtract(u2,u1)
dy        = np.subtract(y2,y1)
Cf        = np.multiply(mu2,np.divide(du,dy))/d_p
Pf        = np.divide(p2,p_i)

with open("Cf_crest.dat",'w') as f:
    for i in range(len(x1)):
        f.write(str('%.7f'%x1[i])+' '+str('%.7f'%Cf[i])+' '\
               +str('%.7f'%Pf[i])+'\n')
    print("Cf data has been output.")
'''


