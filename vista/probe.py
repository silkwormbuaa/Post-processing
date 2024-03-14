# -*- coding: utf-8 -*-
'''
@File    :   probe.py
@Time    :   2024/01/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''
import re
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt
from   .psd              import pre_multi_psd

class Probe:
    
    def __init__( self ) -> None:
        
        """
        fmt:  format of probe data
        rate: sampling rate (how many steps per collection)
        location: probe location
        """
        self.fmt  = 'POINT'
        self.rate = 5 
        self.xyz  = [0.0, 0.0, 0.0]
        self.mode = 'ALL'


# ----------------------------------------------------------------------
# >>> assign a probe                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def assign( self, fmt, rate, xyz, mode ):

        """
        fmt:  format of probe data
        rate: sampling rate (how many steps per collection)
        location: probe location [x,y,z]
        """

        self.fmt  = fmt
        self.rate = rate
        self.xyz  = xyz
        self.mode = mode
        

# ----------------------------------------------------------------------
# >>> write a probe                                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def write( self, file ):
        
        string = f'{self.fmt} {self.rate:3d}  {self.xyz[0]:>13.6e}  ' + \
                 f'{self.xyz[1]:>13.6e}  {self.xyz[2]:>13.6e}  {self.mode}'
        
        file.write( string )
        
       
# ----------------------------------------------------------------------
# >>> CLASS: ProbeFile                                             
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

class ProbeFile:
        
    def __init__( self, filename=None ) -> None:
        
        
        """
        filename : inca_probes.inp file
        
        return:
        self.probes: list of probes
        """
        self.header = '#(1) FORMAT (2) RATE (3) X (4) Y (5) Z (6) MODE'
        self.probes = []
        
        if filename:
            self.read( filename )                                       


# ----------------------------------------------------------------------
# >>> Read a probe file                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read( self, filename ):
        
        with open( filename, 'r' ) as f:
            
            lines = f.readlines()
            
            for line in lines:
                
                if line[0] == '#':
                    
                    continue
                
                else:
                    
                    words = line.split()
                    
                    probe = Probe()
                    
                    probe.assign( words[0], 
                                  int(words[1]), 
                                  [float(words[2]),
                                   float(words[3]),
                                   float(words[4])], 
                                  words[5] )
                    
                    self.probes.append( probe )


# ----------------------------------------------------------------------
# >>> Write a probe file                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def write( self, filename ):
        
        with open( filename, 'w' ) as f:
            
            f.write( self.header )
            
            for probe in self.probes:
                
                f.write( '\n' )
                
                probe.write( f )



# ----------------------------------------------------------------------
# >>> Show probes                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def show( self ):
        
        # plot all probes points in 3D view 
        
        fig = plt.figure()
        ax = fig.add_subplot( projection='3d' )
        
        for probe in self.probes:
            
            ax.scatter( probe.xyz[0], probe.xyz[1], probe.xyz[2] )
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # set aspect equal
        plt.gca().set_aspect('auto', adjustable='box')
        
        plt.show()


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

class ProbeData:
    
    def __init__( self, filename=None, withT=False ):
        
        """
        filename: probe data file
        """
        self.var_list = [ 'step',
                          'time',
                          'u',
                          'v',
                          'w',
                          'rho',
                          'rhoE',
                          'p'    ]
        if withT: self.var_list.append('T')
        
        # probe data dataframe
        self.df = pd.DataFrame()
        
        # psd results dataframe
        self.psd_df = pd.DataFrame()
        
        if filename: self.read( filename)


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read( self, filename ):
        
        with open( filename, 'r' ) as f:
            
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
    

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def cleandata( self, t_start ):
        
        """
        t_start: start time to keep the data
        
        drop the transient data at the beginning
        """

        timelist = np.array( self.df['time'] )
        
        for i in range( len(timelist) ):
            if timelist[i] >= t_start:
                cut_index = i ; break
        
        self.df = self.df[cut_index:]
        

# ----------------------------------------------------------------------
# >>> get fluctuation                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def get_fluc( self, vars ):
        
        """
        vars: list of variables to compute fluctuation
        
        Possible choices ['u','v','w','rho','rhoE','p','T']
        """
        
        for var in vars:
            
            if var not in self.var_list:
                raise ValueError(f"Variable {var} is not in the probe data list.")
            
            mean = np.array( self.df[var] ).mean()
            
            self.df[f'{var}_fluc'] = self.df[var] - mean


# ----------------------------------------------------------------------
# >>> sampling frequency                                         (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

    @property
    def fs( self ):
        
        timelist = np.array( self.df['time'] )
        timespan = timelist[-1] - timelist[0] 
        
        return 1.0/( round( timespan/(len(timelist)-1), 7) )
    
    
# ----------------------------------------------------------------------
# >>> compute probe psd                                   (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def pre_multi_psd( self, var, n_seg, overlap, nfft=None ):

        """
        var: variable to compute PSD
        n_seg: number of segments
        overlap: overlap ratio
        nfft: int, length of the FFT used, if a zero padded FFT is desired. 
              If None, the FFT length is nperseg. Defaults to None.
        
        """
        
        data = np.array( self.df[var] )
        
        freq, pm_psd = pre_multi_psd( data, self.fs, n_seg, overlap, nfft )
        
        if self.psd_df.empty:
            self.psd_df['freq'] = freq
        
        self.psd_df[f'pmpsd_{var}'] = pm_psd



# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

def WriteProbe():

    fname = '/home/wencanwu/my_simulation/STBLI_mid_Re/231124/inca_probes.inp'
    outfile = '/home/wencanwu/my_simulation/STBLI_low_Re/240210/test.inp'
    probes = ProbeFile()
#    probes.read( fname )
#    probes.show()

    xs = np.arange( -30.0, 115, 0.43382353 )
    y  = 0.00001
    z  = 0.00001
    
    for x in xs:
        
        probe = Probe()
        probe.assign( 'POINT', 5, [x,y,z], 'ALL' )
        probes.probes.append( probe ) 
    
    z  = 0.65
    y  = -1.04

    for x in xs:
        
        probe = Probe()
        probe.assign( 'POINT', 5, [x,y,z], 'ALL' )
        probes.probes.append( probe ) 
    
    probes.write( outfile )


# ----------------------------------------------------------------------
# >>> Testing                                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():
    
    filename = '/home/wencanwu/my_simulation/temp/220927_lowRe/probes/probe_x/probe_00097.dat'

    probe = ProbeData( filename, withT=False )
    
    probe.cleandata(20.0)
    
    probe.get_fluc(['p'])
    
    probe.pre_multi_psd('p_fluc', 8, 0.5)
    
    print(probe.psd_df)
    
    
# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()