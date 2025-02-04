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
import pickle
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt
from   vista.timer       import timer
from   vista.colors      import colors as col
from   .psd              import pre_multi_psd
from   .psd              import psd_hann
from   .tools            import define_wall_shape

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
    
    def __init__( self, filename=None, withT=False, step=1 ):
        
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
        
        # probe location
        self.xyz = [0.0, 0.0, 0.0]
        
        # probe index
        self.probe_index = 0
        
        if filename: self.read( filename, step=step )


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

    def read( self, filename, step=1 ):
        
        """
        filename: probe data file, [-9,-4] need to be probe index
        read probe data file and store data into self.df
        """
        
        self.probe_index = int( filename[-9:-4] )
        
        with open( filename, 'r' ) as f:
            
            lines = f.readlines()
            
            # regular expression to read probe location
            
            try:
                x = float( re.search(r'x =(.*?),',lines[0]).group(1) )
                y = float( re.search(r'y =(.*?),',lines[0]).group(1) )
                z = float( re.search(r'z =(.*?)(?:\n)',lines[0]).group(1) )
                self.xyz = [ x, y, z ]
            
            except:
                raise ValueError(f"Probe {self.probe_index:05d}: Cannot find probe location.")
            
            # read in the data body
            
            row = None
            
            for i in range( 1, len(lines), step):
                
                cleanl = lines[i].strip().split()
                
                # skip headers, if there are.
                if cleanl[0].startswith('#'): continue
                
                cleanl = [ parse_float(item,i,self.probe_index) for item in cleanl ]

                if row is None: 
                    
                    row = list()
                    row.append(cleanl)
                
                else:
                    
                    row.append(cleanl)            
  
            self.df = pd.DataFrame( data=row, columns=self.var_list )
    
            # add T if withT is False
            
            if 'T' not in self.var_list:
                self.add_data( 'T', self.df['p']/(self.df['rho']*287.0508571) )

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

    def cleandata( self, t_start, t_end=None ):
        
        """
        t_start: start time to keep the data
        t_end  : end time to keep the data
        
        drop the transient data at the beginning or keep the data between start 
        and end.
        """

        timelist = np.array( self.df['time'] )
        
        if t_start > timelist[-1]:
            raise ValueError(f"Probe {self.probe_index:05d}: Start time {t_start} is out of probe data range.")
        
        for i in range( len(timelist) ):
            if timelist[i] >= t_start:
                cut_index = i ; break
        
        if t_end is not None:
            
            for i in range( cut_index, len(timelist) ):
                if timelist[i] >= t_end: break; 
            end_index=i
            self.df = self.df[cut_index:end_index]
            
        else:
            self.df = self.df[cut_index:]
        

# ----------------------------------------------------------------------
# >>> add_data                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/09/18  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def add_data( self, varname, dataarray ):

        """
        varname: name of the variable
        dataarray: array of data
        
        add data to the dataframe
        """
        
        self.var_list.append( varname )
        self.df[varname] = dataarray
        

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
# >>> compute probe pre-multiplied psd                                   (Nr.)
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

    def pre_multi_psd( self, var, n_seg, overlap, nfft=None, normalize=True ):

        """
        var: variable to compute PSD
        n_seg: number of segments
        overlap: overlap ratio
        nfft: int, length of the FFT used, if a zero padded FFT is desired. 
              If None, the FFT length is nperseg. Defaults to None.
        
        """
        
        data = np.array( self.df[var] )
        
        freq, pm_psd = pre_multi_psd( data, self.fs, n_seg, overlap, nfft, normalize )
        
        if self.psd_df.empty:
            self.psd_df['freq'] = freq
        
        self.psd_df[f'pmpsd_{var}'] = pm_psd


# ----------------------------------------------------------------------
# >>> compute probe psd                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/04/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_psd( self, var, n_seg, overlap, nfft=None ):
        
        """
        var: variable to compute PSD
        n_seg: number of segments
        overlap: overlap ratio
        nfft: int, length of the FFT used, if a zero padded FFT is desired. 
              If None, the FFT length is nperseg. Defaults to None.
        
        """
        
        data = np.array( self.df[var] )
        
        freq, psd = psd_hann( data, self.fs, n_seg, overlap, nfft )
        
        if self.psd_df.empty:
            self.psd_df['freq'] = freq
        
        self.psd_df[f'psd_{var}'] = psd 

# ----------------------------------------------------------------------
# >>> write psd dataframe                                        (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/15  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def write_psd( self, filename=None ):
        
        """
        dump self.xyz and seld.psd_df
        """
        
        if filename is None:
            filename = f'psd_{self.probe_index:05d}.pkl'
        
        with open( filename, 'wb' ) as f:
            
            pickle.dump( self.xyz, f)
            pickle.dump( self.psd_df, f )


# ----------------------------------------------------------------------
# >>> time_index                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/09/15  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def time_index( self, ts: np.array, buff=0.001 ) -> np.array:
        
        """
        ts: array of CONSECUTIVE time points\n
        
        return: array of index corresponding to nearest time points
        """
        
# ----- check if the give time range is within the probe data range

        t_min = self.df['time'].min()
        t_max = self.df['time'].max()
        
        if (t_min-buff) > ts[0] or (t_max+buff) < ts[1]:
            raise ValueError(f"Time range {ts} is out of probe data covered range [{t_min},{t_max}].")

# ----- prepare data
        
        timep = np.array( self.df['time'] )
        index = list()
        i = 0
        
# ----- loop over every snapshot time point

        for t in ts:

            # if (timep[i] < t) is put forward, the index will overfloat
            # so we need to check if i is out of range first.
            while i < len(timep) and timep[i] < t : i += 1
            
            # probe data may not cover the last time point, but the buffer 
            # already ensures the last time point is within the probe data range
            # so we just need to add the last index.
            if   i == len(timep): 
                index.append(i-1)
                
            elif timep[i] >= t:    
                index.append(i)
                
        return np.array(index)
        

# ----------------------------------------------------------------------
# >>> read psd                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/15  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read_psd( self, filename ):
        
        """
        load self.xyz and seld.psd_df
        """
        
        self.probe_index = int( filename[-9:-4] )
        
        with open( filename, 'rb' ) as f:
            
            self.xyz = pickle.load( f )
            self.psd_df = pickle.load( f )


# ----------------------------------------------------------------------
# >>> parse float considering bad values                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/15  - created
#
# Desc
#    some probe data contains bad values, give a warning to user if so.
# ----------------------------------------------------------------------

def parse_float(value, i, index):
    
    try:
        return float(value)
    
    except ValueError:
        print(col.fg.yellow,f"Probe {index:05d}: Bad value {value} is found in row {i}.",col.reset)
        return 0.0



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
    
    filename = '/home/wencanwu/temp/probe/probe_00001.dat'

    with timer("Read probe data"):
        probe = ProbeData( filename, withT=True )
    
    probe.cleandata(15.0)

"""
# test finding index from given time points

    ts = np.linspace(20.0, 61.0, 4101)
    index = probe.time_index( ts )
    
    df = pd.DataFrame({'ts':ts, 'index':index})
    
    df.to_csv('/home/wencanwu/test/probe_find_index/index.csv', index=False)
    
    # extract the data at given time points
    
    df_extract = probe.df.iloc[index]
    df_extract.reset_index(drop=True,inplace=True)
    print(df_extract)

"""    
    
#    probe.get_fluc(['p'])
    
#    probe.pre_multi_psd('p_fluc', 8, 0.5)
    
#    print(probe.psd_df)
    
    
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