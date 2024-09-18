# -*- coding: utf-8 -*-
'''
@File    :   species.py
@Time    :   2024/09/18 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import numpy             as     np


# ----------------------------------------------------------------------
# >>> get_visc                                                 (Nr.)
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

def get_visc( ts, law='P07', case='midRe' ):
    
    """
    ts     : static temperature
    
    law    : 'sutherlands', 'P05', or 'P07
    case   : 'lowRe' or 'midRe'
    return : mu
    """
    
    u_ref      = 507.0
    l_ref      = 1.0
    T_ref      = 160.15
    rho_ref    = 0.9886
    T_ref_sl   = 273.15
    mu_ref_sl  = 17.16e-6
    S_ref_sl   = 110.4
    
    if case == 'lowRe':   Re_ref = 2230.0
    elif case == 'midRe': Re_ref = 9636.0
    
    visc_ref = rho_ref * l_ref * u_ref / Re_ref

    def sutherland( t ):
        
        mu = mu_ref_sl * (t/T_ref_sl)**1.5 * (T_ref_sl + S_ref_sl) / (t + S_ref_sl)
        
        return mu
    
    def P05( t ):
        
        p_ref = visc_ref / T_ref**0.5
        mu = p_ref * abs(t)**0.5
        
        return mu
        
    def P07( t ):
        
        p_ref = visc_ref / T_ref**0.7
        mu = p_ref * abs(t)**0.7
        
        return mu

    def process( ts, method=law ):
        
        if method == 'P05':
            return P05(ts)
        elif method == 'P07':
            return P07(ts)
        elif method == 'sutherland':
            return sutherland( ts )
        else:
            raise ValueError("method should be 'P07' or 'sutherland'.")
        

    if isinstance( ts, float ) or isinstance( ts, int):
        return process( ts )
    elif isinstance( ts, np.ndarray ):
        return process( ts)
    elif isinstance( ts, list ):
        return [ process( t ) for t in ts ]
    else:
        raise ValueError("ts should be a float, int, list or np.ndarray.")


# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
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

def Testing():

    t1 = 223.03
    t2 = [250.0,260.0,270.0,280.0]
    t3 = np.array([250.0,260.0,270.0,280.0])
    
    law  = 'P05'
    case = 'lowRe'
    
    print( get_visc( t1, law=law, case=case ) )
    print( get_visc( t2, law=law, case=case ) )
    print( get_visc( t3, law=law, case=case ) )



# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
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

if __name__ == "__main__":

    Testing()
