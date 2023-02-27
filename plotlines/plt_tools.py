#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_tools.py
@Time    :   2022/10/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Functions which will facilitate plotting 
'''
import re

import numpy             as     np

import pandas            as     pd

import matplotlib.pyplot as     plt

# ----------------------------------------------------------------------
# >>> CLASS - DATAFRAME FOR PLOTTING                            ( 1-0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/17  - created
#
# Desc
#
# ----------------------------------------------------------------------


class PlotDataframe():

# ----------------------------------------------------------------------
# >>> Initialization by reading in data                         ( 1-1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/17  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def __init__( self, filename ):
        
        with open( filename, 'r' ) as f:
            
            # read in the header as columns variables
            #   str.strip() remove leading and trailing whitespaces
            #   str.split() split a string into a string list, separator can 
            #               be specified, default is whitespace. 
            
            self.var_list = f.readline().strip().split()
        
            row = None
        
            line = f.readline()
            while line:
                
                cleanl = line.strip().split()
                
                # transfer from str to float
                
                cleanl = [float(item) for item in cleanl]
                
                if row is None:
                    row = np.array(cleanl)
                else: 
                    row = np.vstack((row, cleanl))
                    
                line = f.readline()
            
            self.df = pd.DataFrame( data=row, columns=self.var_list )



# ----------------------------------------------------------------------
# >>> FETCH DATA FOR X,Y AXIS                                  ( 1-2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/17  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def fetch_xy( self, var_x, var_y ):

        value_x = np.array( self.df[var_x] )
        
        value_y = np.array( self.df[var_y] )
        
        return ( value_x, value_y ) 
    

# ----------------------------------------------------------------------
# >>> READ NORMALIZATION VALUES                                ( 1-3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/17  - created
#
# Input
#
# - statistic file containing values computed from IB force
#
# ----------------------------------------------------------------------

    def read_norm( self, filename ):
        
        with open( filename, 'r' ) as f:
            
            # skip header
            
            line = f.readline()
            
            line = f.readline().strip().split()
            self.tau_w = float( line[0] )
            self.rho_w = float( line[1] )
            self.u_tau = float( line[2] )
            self.nu    = float( line[3] )
            self.mu    = float( line[4] )
            self.lv    = float( line[5] )


# ----------------------------------------------------------------------
# >>> NORMALIZE VALUES                                         ( 1-4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/17  - created
#
# Desc
#  
# - normalize profile data(y+,u+,u'v' etc) with data from read_norm()
#
# ----------------------------------------------------------------------

    def get_norm( self ):
        
        # get normalized values
        
        y_plus = np.array(self.df['y']) / self.lv
        
        y_s_plus = np.array(self.df['y_s']) / self.lv

        u_plus   = np.array(self.df['<u>']) / self.u_tau
        
        rho = np.array( self.df['<rho>']  )
        uu  = np.array( self.df['<u`u`>'] )
        vv  = np.array( self.df['<v`v`>'] )
        ww  = np.array( self.df['<w`w`>'] )
        uv  = np.array( self.df['<u`v`>'] )

        uu_plus  = np.multiply(rho,uu) / abs(self.tau_w)
        vv_plus  = np.multiply(rho,vv) / abs(self.tau_w)
        ww_plus  = np.multiply(rho,ww) / abs(self.tau_w)
        uv_plus  = np.multiply(rho,uv) / abs(self.tau_w)
        
        # add normalized results to dataframe
        
        self.df['y_plus']   = y_plus 
        self.df['y_s_plus'] = y_s_plus
        self.df['u_plus']   = u_plus
        
        self.df['<u`u`>+']  = uu_plus
        self.df['<v`v`>+']  = vv_plus
        self.df['<w`w`>+']  = ww_plus
        self.df['<u`v`>+']  = uv_plus

        # normalize temperature profile 
        
        self.T_inf = 160.15
        
        T    = np.array( self.df['<T>'] )
        T_nd = T/self.T_inf        
        
        self.df['T_nd'] = T_nd 
        

# ----------------------------------------------------------------------
# >>> SHIFT Y COORDINATE                                       ( 1-5 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/10/17  - created
#
# Desc
#
# - shift y coordinate based on roughness average elevation
#
# ----------------------------------------------------------------------

    def shift_y( self, dy ):
        
        y = np.array( self.df['y'] )
        
        y_s = np.add( y, dy )

        self.df['y_s'] = y_s
        

# ----------------------------------------------------------------------
# >>> SHIFT X COORDINATE                                       ( 1-6 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/01  - created
#
# Desc
#
# - shift x coordinate from x_imp
#
# ----------------------------------------------------------------------

    def shift_x( self, x_imp, delta=None ):
        
        x = np.array( self.df['x'] )
        
        x_s = np.subtract( x, x_imp )
        
        self.df['x_s'] = x_s 
        
        if delta is not None:
            
            self.df['x_s'] = np.divide( x_s, delta)

  
# ----------------------------------------------------------------------
# >>> van Driest transform                                     ( 1-7 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/12/11  - created
#
# Desc
#
# - van Driest transform mean velocity^+
# - Must do 
#   (1) read_norm( statistic_average.dat ) to have rho_w
#   (1) get_norm() to have u+ first
# ----------------------------------------------------------------------

    def vd_transform( self, mode="bottom" ):
        
        rho = np.array( self.df['<rho>'] )
        u_plus = np.array( self.df['u_plus'])
        
        # if integrate from the crest ( y=0 )
        if mode == "crest" :
            
            y0_index = int( self.df[self.df['y']==0].index[0] ) 
            
            rho[:y0_index-1] = 0
            u_plus[:y0_index] = 0

        # by default, integrate from the bottom ( wall_dist=0 )
        if self.rho_w is None:
            
            raise ValueError("Please read_norm() to have rho_w!")
        
        else:
                        
            rho_ratio = np.sqrt( rho/self.rho_w )
        
        u_plus_vd = np.zeros( np.size(u_plus) )
        
        for i in range( len(u_plus_vd) ):
            
            u_plus_vd[i] = np.trapz( rho_ratio[0:i+1],u_plus[0:i+1] )
        
        self.df['u_plus_vd'] = u_plus_vd      
            
# ----------------------------------------------------------------------
# >>> Get turbulent Mach number                                 ( 1-8 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/01/29  - created
#
# Desc
#
# - get the turbulent Mach number data
#
# ----------------------------------------------------------------------
    def get_Mt( self ):
        
        uu = np.array( self.df['<u`u`>'] )
        vv = np.array( self.df['<v`v`>'] )
        ww = np.array( self.df['<w`w`>'] )
        
        #Vp - fluctuating velocity
        Vp = np.sqrt( uu + vv + ww )
        
        T = np.array( self.df['<T>'] )
        gamma = 1.4
        R = 287.06
        c = np.sqrt(gamma*R*T)
        
        # avoid c = 0 as dominator
        Mt = np.zeros( (np.size(c)) )
        for i in range( np.size(c)) :
            if c[i] > 0.0:
                Mt[i] = Vp[i]/c[i]
        
        self.df['Mt'] = Mt
        
           
# ----------------------------------------------------------------------
# >>> Read psd data                                            ( 2-1 )
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
# - reading psd_00xxx.dat file
#
# ----------------------------------------------------------------------

def read_psd( filename ):
    
    with open( filename, 'r' ) as f:
        
        lines = f.readlines()
        
        # regular expression to read probe location
        
        xp = float( re.search(r'x=(.*?) y',lines[0]).group(1) )
        yp = float( re.search(r'y=(.*?) z',lines[0]).group(1) )
        zp = float( re.search(r'z=(.*?)(?:\n)',lines[0]).group(1) )
        
        # var_list 
        
        var_list = [
            'freq',
            'pprime_fwpsd',
            'St',
            'nd_pprime_fwpsd'
        ]
        
        # read in data body
        
        row = None
        
        for i in range( 1, len(lines) ):
            
            cleanl = lines[i].strip().split()
            
            cleanl = [ float(item) for item in cleanl ]
            
            if row is None:
                
                row = list()
                row.append( cleanl )
                
            else:
                
                row.append( cleanl )
        
    df = pd.DataFrame( data=row, columns=var_list )
     
    St = np.array( df['St'] ).tolist()
    nd_fwpsd = np.array( df['nd_pprime_fwpsd'] ).tolist()
    x = (np.array([xp] * len(St))-50.4)/5.2  
    # get a xp array,length equals the number of frequency
    
    return ( x, St, nd_fwpsd )
        
        
        