# -*- coding: utf-8 -*-
'''
@File    :   line.py
@Time    :   2023/10/03 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   LineData, mainly for easing plotting lines
'''

import numpy             as     np
import pandas            as     pd
from   scipy             import fft

from   .tools            import find_indices


class LineData:
# ----------------------------------------------------------------------
# >>> LineData                                               (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def __init__( self, filename=None, df=None ):
        
        """
        filename : file containing header and data
        
        if given filename, initialize from reading file; 
        else generate an empty LineData object.
        """
        
        self.df = None
        
        self.label = None
        
        self.color = None
        
        self.width = None
        
        self.lstyle = None
        
        if filename is not None:
            
            self.read( filename )
        
        if df is not None:
            
            self.df = df
            self.vars = df.columns
    
    
# ----------------------------------------------------------------------
# >>> Read from file                                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def read( self, filename ):
        
        """
        filename : file containing header and data
        
        return : self.df 
        """
        
        with open( filename, 'r' ) as f:
            
            self.vars = f.readline().strip().split()
            
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
            
            self.df = pd.DataFrame( data=row, columns=self.vars )


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/09  - created
#
# Desc
#
# ----------------------------------------------------------------------
    
    def fft_kz( self ):
        
        """
        x_str  : 'x','y','z'
        y_str : 'u`','v`','w`' or other variable name to be fft-ed
        
        return : wavenumber k and y_k
        """
        
        z = np.array( self.df['z'] )
        u_fluc = np.array( self.df['u`'] )
        v_fluc = np.array( self.df['v`'] )
        w_fluc = np.array( self.df['w`'] )
        E_fluc = u_fluc**2 + v_fluc**2 + w_fluc**2

        # 1-D fourier transform
        
        u_k = fft.fft( u_fluc )
        v_k = fft.fft( v_fluc )
        w_k = fft.fft( w_fluc )
        E_k = fft.fft( E_fluc )
        
        # wave number 
        
        k_z = 2.0 * np.pi * fft.fftfreq( len(z), z[1]-z[0]) 

        self.df['k_z'] = k_z
#        self.df['u_k'] = u_k
#        self.df['v_k'] = v_k 
#        self.df['w_k'] = w_k 
 
        self.df['E_uu'] = np.abs(u_k)*2
        self.df['E_vv'] = np.abs(v_k)*2
        self.df['E_ww'] = np.abs(w_k)*2
        self.df['E_k']  = np.abs(E_k)
        

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/25  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def sparse_log( self, col, spacing, base=10):
        
        """
        df: input pd.DataFrame
        col: col name
        spacing: spacing for sparsed results
        base: base of log, default is 10.
        
        return: sparsed dataframe
        """

        index = []
        
        x = np.array( self.df[col] )
        
        # check if x is monotonic increasing

        if np.all( np.diff(x) > 0 ):
            
            # find the first element of x
            
            index += [0,1]
            
            # find the rest elements of x
            
            for i in range(2,len(x)):
                if np.emath.logn( base, x[i]/x[index[-1]] ) >= spacing:
                    index.append(i)
        
        
        return self.df.iloc[index]


class ProfileData( LineData ):
# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/03  - created
#
# Desc
#
# ----------------------------------------------------------------------
    
    def __init__( self, filename=None ):
        
        super().__init__( filename )
    
        self.tau_ave = None
        self.rho_ave = None
        self.mu_ave  = None
        self.u_tau   = None
        self.lv      = None


# ----------------------------------------------------------------------
# >>> shift y                                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def shift_y( self, dy ):
        
        y = np.array( self.df['y'] )
        
        self.df['ys'] = y + dy


# ----------------------------------------------------------------------
# >>> Read norm                                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/03  - created
#
# Desc
#   - read wall statistics from file
#
# ----------------------------------------------------------------------

    def inner_normalize( self, filename ):
        
        """
        filename : 'wall_statistics.dat' file
        
        return : y+,ys+,u+,u`u`+(density scaled),...,u`v`+ in self.df
        """
        
        with open( filename, 'r' ) as f:
            
            lines = f.readlines()
            
            self.tau_ave = float( lines[0].strip().split()[1] )
            self.rho_ave = float( lines[1].strip().split()[1] )
            self.mu_ave  = float( lines[2].strip().split()[1] )
            self.u_tau   = float( lines[3].strip().split()[1] )
            self.lv      = float( lines[4].strip().split()[1] )
    
        self.df['y+']  = np.array( self.df['y'] ) / self.lv
        self.df['ys+'] = np.array( self.df['ys'] ) / self.lv
        self.df['u+']  = np.array( self.df['u']) / self.u_tau
        
        rho = np.array( self.df['rho'] )
        
        self.df['u`u`+'] = rho*np.array( self.df['u`u`'] ) / self.tau_ave
        self.df['v`v`+'] = rho*np.array( self.df['v`v`'] ) / self.tau_ave
        self.df['w`w`+'] = rho*np.array( self.df['w`w`'] ) / self.tau_ave
        self.df['u`v`+'] = rho*np.array( self.df['u`v`'] ) / self.tau_ave


# ----------------------------------------------------------------------
# >>> van Driest transform                                          (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def vd_transform( self, mode="bottom" ):
        
        """
        mode : by default 'bottom', 'crest' is wrong
        return : self.df['u+_vd']
        
        """
        
        rho = np.array( self.df['rho'] )
        u_plus = np.array( self.df['u+'])
        
        # if integrate from the crest ( y=0 )
        if mode == "crest" :
            
            y0_index = int( self.df[self.df['y']==0].index[0] ) 
            
            rho[:y0_index-1] = 0
            u_plus[:y0_index] = 0

        # by default, integrate from the bottom ( wall_dist=0 )
        if self.rho_ave is None:
            raise ValueError("Please inner_normalize() to have rho_ave!")
        else:
            rho_ratio = np.sqrt( rho/self.rho_ave )
        
        u_plus_vd = np.zeros( np.size(u_plus) )
        
        for i in range( len(u_plus_vd) ):
            u_plus_vd[i] = np.trapz( rho_ratio[0:i+1],u_plus[0:i+1] )
            
        
        self.df['u+_vd'] = u_plus_vd    


# ----------------------------------------------------------------------
# >>> Compute displacement thickness                             (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/28  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def compute_boundary_layer_property( self, u_ref ):
        
        """
        u_ref : reference velocity
        
        return : delta_99, delta_star(displacement), theta(momentum)
        """
        
        u = np.array(self.df['u'])
        y = np.array(self.df['y'])
        rho = np.array(self.df['rho'])
        
        il,ir = find_indices( u, u_ref*0.99, mode='sequential' )
        
        rho_e = rho[ir]
        
        # linear interpolate between left and right points
        
        delta = (u_ref*0.99 - u[il]) / (u[ir] - u[il]) * (y[ir] - y[il]) + y[il]
        
        rho_e = (rho[ir]-rho[il]) / (y[ir]-y[il]) * (delta - y[il]) + rho[il]
        
        u[ir] = u_ref*0.99
        y[ir] = delta
        rho[ir] = rho_e
        
        # compute displacement thickness and momentum thickness
        
        delta_star = np.trapz( ( 1.0 - u[:ir]*rho[:ir]/u_ref/rho_e ), y[:ir] )
        
        theta = np.trapz( ( 1.0 - u[:ir]/u_ref ) * rho[:ir]*u[ir]/rho_e/u_ref,
                          y[:ir] )
        
        return delta, delta_star, theta
        
  
        
# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    file = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results/profile/profile_mean.dat'
    normfile = '/media/wencanwu/Seagate Expansion Drive/temp/221221/results/profile/wall_statistics.dat'
    
    line1 = ProfileData( file )
    
    print(line1.df)
    
    line1.shift_y( 0.26 )

    line1.inner_normalize( normfile )
    
    line1.vd_transform()
    
    print(line1.df)


# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/10/03  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
