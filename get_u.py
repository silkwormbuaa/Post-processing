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
FoldPath = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/TP_stat"
OutPath  = "/home/wencanwu/my_simulation/temp/Low_Re_Luis/DataPost/"

zonegrp = GetZonegrp(FoldPath)
with timer("get lines"):
    line_loc = [-60.0, 0.0, -60.0, 10.0]    
    GetLine(line_loc,zonegrp,FoldPath,OutPath)
    line_loc = [-40.0, 0.0, -40.0, 10.0]    
    GetLine(line_loc,zonegrp,FoldPath,OutPath)


#        print(dataset.num_solution_times)
#        print(dataset.VariablesNamedTuple)
#        print(dataset.num_variables)
#        print(dataset.num_zones)
#        print(dataset.title)
#        print(dataset.variable_names)
