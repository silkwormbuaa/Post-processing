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

import numpy             as     np

import pandas            as     pd

import matplotlib.pyplot as     plt

# ----------------------------------------------------------------------
# >>> CLASS - DATAFRAME FOR PLOTTING                               ( 0 )
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
# >>> Initialization by reading in data                           ( 1 )
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
# >>> FETCH DATA FOR X,Y AXIS                                     ( 2 )
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

        value_x = np.array(self.df[var_x])
        
        value_y = np.array(self.df[var_y])
        
        return (value_x,value_y) 
    

# ----------------------------------------------------------------------
# >>> READ NORMALIZATION VALUES                                   ( 3 )
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

    def read_norm(self, filename):
        
        with open(filename, 'r') as f:
            
            # skip header
            
            line = f.readline()
            
            line = f.readline().strip().split()
            self.tau_w = float(line[0])
            self.rho_w = float(line[1])
            self.u_tau = float(line[2])
            self.nu    = float(line[3])
            self.mu    = float(line[4])
            self.lv    = float(line[5])

# ----------------------------------------------------------------------
# >>> SHIFT Y COORDINATE                                          ( 4 )
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

    def shift_y( self, dy ):
        
        y = np.array(self.df['y'])
        
        y_s = np.add( y, dy )

        self.df['y_s'] = y_s


# ----------------------------------------------------------------------
# >>> NORMALIZE VALUES                                            ( 5 )
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

    def get_norm( self ):
        
        # get normalized values
        
        y_plus = np.array(self.df['y']) / self.lv
        
        y_s_plus = np.array(self.df['y_s']) / self.lv

        u_plus   = np.array(self.df['<u>']) / self.u_tau
        
        rho = np.array(self.df['<rho>'])
        uu = np.array(self.df['<u`u`>'])
        vv = np.array(self.df['<v`v`>'])
        ww = np.array(self.df['<w`w`>'])
        uv = np.array(self.df['<u`v`>'])
        
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

    
'''  #test scripts  
file = '/home/wencanwu/my_simulation/temp/220927_lowRe/mean_result_ib.dat' 
norm_file = '/home/wencanwu/my_simulation/temp/220927_lowRe/statistic_average.dat'

u_1 = PlotDataframe(file)

print(u_1.df)

x,y = u_1.fetch_xy('y','<u>')
   
u_1.shift_y(1)

u_1.read_norm(norm_file)

u_1.get_norm()



print(u_1.df)

'''
#y,u = fetch_xy(df, 'y', '<u>')

#print(x)
#print(y)



#fig, ax = plt.subplots(figsize=[10,8])
#
#ax.plot(u_1.df['y'],u_1.df['<u>'])
#
#plt.show()
#print(type(y[0]))