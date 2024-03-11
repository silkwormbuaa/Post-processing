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

import matplotlib.pyplot as     plt


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
    
    def __init__( self, datafile, withT=True ):
        
        """
        datafile: probe data file
        withT: new INCA since Jan 2023, set True by default. Else set False.
        """
        
        self.probe_index = int( datafile[-9:-4] )
        
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
        if withT:                          # updated inca since Jan 2023,
            self.var_list.append('T')      # by default, T is added.     
                                             

        with open( datafile, 'r' ) as f:

            lines = f.readlines()
            
            # regular expression to read probe location
            
            self.x = float( re.search(r'x =(.*?),',lines[0]).group(1) )
            self.y = float( re.search(r'y =(.*?),',lines[0]).group(1) )
            self.z = float( re.search(r'z =(.*?)(?:\n)',lines[0]).group(1) )
            
            # read in the data body
            
            row = None
            
            for i in range( 1, len(lines) ):
                
                cleanl = lines[i].strip().split()
                
                cleanl = [ float(item) for item in cleanl ]

                if row is None: 
                    
                    row = list()
                    row.append(cleanl)
                
                else:
                    
                    row.append(cleanl)            
  
            self.df = pd.DataFrame( data=row, columns=self.var_list )

        # set case related constant parameters
        # - boundary layer thickness; freestream velocity
        
        self.delta_0 = 5.2
        self.U_inf   = 507
#        self.Lsep    = 9.805522
        
        
        
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
# ----------------------------------------------------------------------

    def cleandata( self, starttime ):

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
        
        self.freq_sample = 1.0/self.meaninterval
        
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
            pprime,                fs=self.freq_sample,  
            window='hann',         nperseg=self.nperseg, 
            noverlap=self.noverlap,scaling='density',
            nfft=len(pprime))
        
        self.mean_power = np.sum(self.pprime_psd * self.freq[1])
        
        self.pprime_fwpsd = np.multiply(self.freq,self.pprime_psd) 
        
        self.nd_pprime_fwpsd = self.pprime_fwpsd / self.mean_power
        
        self.St = self.freq * self.delta_0 / self.U_inf
        
#        self.St_Lsep = self.freq * self.Lsep / self.U_inf
        
 # ----------------------------------------------------------------------
 # >>> Plot PSD                                                ( 3 )
 # ----------------------------------------------------------------------
 #
 # Wencan Wu : w.wu-3@tudelft.nl
 #
 # History
 #
 # 2022/11/13  - created
 #
 # Desc
 #
 # ----------------------------------------------------------------------
        
    def plot_psd( self, Lsep ):
        
        fig, ax = plt.subplots( figsize=[10,8], constrained_layout=True )
        
        St_Lsep = self.St * Lsep
        
#        ax.plot( self.St,
#                 self.nd_pprime_fwpsd )
        
        ax.minorticks_on()
        
#        ax.set_xscale( 'symlog', linthresh = 1)

        ax.semilogx(St_Lsep, self.nd_pprime_fwpsd,'k', linewidth=1)
        
        ax.set_xlim( [0.01,100] )
        
#        x_minor = matplotlib.ticker.LogLocator( 
#                            base = 10.0, subs = np.arange(0.001,10.0) )
#        
#        ax.xaxis.set_minor_locator( x_minor )
        
        plt.show()
        

# ----------------------------------------------------------------------
# >>> Write PSD                                                ( 4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/14  - created
#
# Desc
#
# - write PSD data into a file 
# - with header of location
#
# ----------------------------------------------------------------------

    def write_psd( self ):
        
        outfile = 'psd_%05d.dat'%self.probe_index
        
        with open( outfile, 'w') as f:
            
            f.write( 'x=%15.6e  y=%15.6e  z=%15.6e'%(self.x,self.y,self.z) )
            f.write('\n')

            for i in range(len(self.freq)):
                f.write(str('{:<17.8e}'.format(self.freq[i])))
                f.write(str('{:<17.8e}'.format(self.pprime_fwpsd[i])))
                f.write(str('{:<17.8e}'.format(self.St[i])))
                f.write(str('{:<17.8e}'.format(self.nd_pprime_fwpsd[i])))
                f.write('\n')
        
        
        
        
        
        
