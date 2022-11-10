#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_PSD.py
@Time    :   2022/11/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   For plotting power spectrum density
'''

import os 

import re

import numpy             as     np

import pandas            as     pd


# ----------------------------------------------------------------------
# >>> Class - Probe Data                                           ( 0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/10  - created
#
# Desc
#
# - One probe data file corresponds to one instance
#
# ----------------------------------------------------------------------

class ProbeData:
    
    def __init__( self, datafile ):
        
        self.var_list = [
                          'step',
                          'time',
                          'u',
                          'v',
                          'w',
                          'rho',
                          'rhoE',
                          'p'                          
                        ]

        with open( datafile, 'r' ) as f:
            
            firstline = f.readline()
            
            # regular expression to read probe location
            
            self.x = float( re.search(r'x =(.*?),',firstline).group(1) )
            self.y = float( re.search(r'y =(.*?),',firstline).group(1) )
            self.z = float( re.search(r'z =(.*?)(?:\n)',firstline).group(1) )

            row = None
            
            line = f.readline()

            while line:
                
                cleanl = line.strip().split()
                
                cleanl = [float(item) for item in cleanl]
                
                if row is None: 
                    row = np.array(cleanl)
                else:
                    row = np.vstack((row,cleanl))
                
                line = f.readline()

                if int(cleanl[0])%1000 == 0:
                    print("current step:%s",cleanl[0])
                
            self.df = pd.DataFrame( data=row, columns=self.var_list )
            
            