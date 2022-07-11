# -*- coding: utf-8 -*-
"""
Original code credit to Weibo
"""
from pyexpat.errors import XML_ERROR_INCOMPLETE_PE
import tecplot as tp
#from tecplot.constant import *
import pandas as pd
import os
from os import path
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
from ReadIn import *
from timer import timer
#from time import time
import logging as log
from BLprofile import *

log.basicConfig(level=log.INFO)

#%% Read plt data from INCA
FoldPath = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/TP_stat"
OutPath  = "/home/wencanwu/my_simulation/temp/090522_lowRe_256/DataPost/"
zonegrp2 = GetZonegrp(FoldPath)
with timer("get lines"):
    line_loc = [-58.25, -0.52, -58.25, 30.0]    
    GetLine(line_loc,zonegrp2,FoldPath,OutPath,1)
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
"""
FoldPath = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/TP_stat"
OutPath  = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/DataPost/"

zonegrp = GetZonegrp(FoldPath)
with timer("get lines"):
    line_loc = [-58.25, 0.0, -58.25, 30.0]    
    GetLine(line_loc,zonegrp,FoldPath,OutPath,2)
"""