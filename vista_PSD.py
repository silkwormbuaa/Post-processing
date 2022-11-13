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

import re

import numpy             as     np

import pandas            as     pd

from   scipy             import signal


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

            lines = f.readlines()
            
            # regular expression to read probe location
            
            self.x = float( re.search(r'x =(.*?),',lines[0]).group(1) )
            self.y = float( re.search(r'y =(.*?),',lines[0]).group(1) )
            self.z = float( re.search(r'z =(.*?)(?:\n)',lines[0]).group(1) )
            
            # read in the data body
            
            row = None
            
            for i in range(1,len(lines)):
                
                cleanl = lines[i].strip().split()
                
                cleanl = [float(item) for item in cleanl]

                if row is None: 
                    
                    row = list()
                    row.append(cleanl)
                
                else:
                    
                    row.append(cleanl)            
  
            self.df = pd.DataFrame( data=row, columns=self.var_list )


# ----------------------------------------------------------------------
# >>>  Clean Data                                              ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/11  - created
#
# Desc
#
# - drop transient data at the beginning
# - trim data to multiples of N of segments points for PSD(welch)
# - get mean time interval, frequency, mean p and p'
#
#
# ----------------------------------------------------------------------

    def cleandata( self, starttime, Nseg=None):

#        - this line can also work, but slower
#        self.df = self.df.drop(self.df[self.df['time']<starttime].index)

#       drop transient data before starttime

        timelist = np.array(self.df['time'])
        
        for i in range(len(timelist)):
            
            if timelist[i] > starttime:
                
                cut_index = i
                break

        self.df = self.df[cut_index:]
        
#       drop data to meet Welch's method(has to be multiple of segments)        
#
#        nr_drop = len(self.df)%Nseg
#        
#        print("now drop last %d lines data to meet Welch's method."%nr_drop)
#        
#        self.df = self.df[0:-nr_drop]
#              
#        self.df = self.df.reset_index()
        
#       get mean time interval and frequency

        timespan = timelist[-1] - timelist[0] 
        
        self.meaninterval = round( timespan/(len(timelist)-1), 7 )
        
        self.frequency = 1.0/self.meaninterval
        
#       get mean pressure and p'

        p = np.array( self.df['p'] )

        self.meanp = np.mean( p ) 
        
        pprime = np.subtract( p, self.meanp )
        
        self.df['pprime'] = pprime


# ----------------------------------------------------------------------
# >>> Welch Method                                              ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/11  - created
#
# Desc
#
# - use scipy library to apply Welch Method
#
# ----------------------------------------------------------------------

    def psd( self, n_seg, p_overlap ):
        
        pprime = np.array( self.df['pprime'] )
        
        self.nperseg = int( len(pprime)//(1+(n_seg-1)*(1-p_overlap)) )
        
        self.noverlap = int( self.nperseg*p_overlap )
                
        self.freq, self.pprime_psd = signal.welch(
            pprime,       fs=self.frequency,  window='hann', 
            nperseg=self.nperseg, noverlap=self.noverlap )
        
        
        
        pass
        
        
        
        
        
