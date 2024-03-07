# -*- coding: utf-8 -*-
'''
@File    :   probe.py
@Time    :   2024/01/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import numpy             as     np
import matplotlib.pyplot as     plt


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

def Testing():

    fname = '/home/wencanwu/my_simulation/STBLI_mid_Re/231124/inca_probes.inp'
    outfile = '/home/wencanwu/my_simulation/STBLI_low_Re/240210/test.inp'
    probes = ProbeFile()
#    probes.read( fname )
#    probes.show()

    xs = np.arange( -20.0, 115, 0.43382353 )
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