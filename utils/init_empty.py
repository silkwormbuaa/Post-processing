# -*- coding: utf-8 -*-
'''
@File    :   init_empty.py
@Time    :   2023/05/08 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import numpy             as np

# ----------------------------------------------------------------------
# >>> Initialize empty float array of given dimensions           ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/08  - adapted from L.Laguarda
#
# Desc
#  
# ----------------------------------------------------------------------

def init_1Dflt_empty( elem, precision ):

   if   precision == 4:
      return np.empty( elem, dtype='f4' )

   else:
      return np.empty( elem, dtype='f8' )



# ----------------------------------------------------------------------
# >>> Initialize empty float 2D-array of given                   ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/08  - adapted from L.Laguarda
#
# Desc
#  - row-majored : order='C'
# ----------------------------------------------------------------------

def init_2Dflt_empty( rows, cols, precision ):

   if precision == 4:
      return np.empty( (rows,cols), dtype='f4', order='C' )

   else:
      return np.empty( (rows,cols), dtype='f8', order='C' )
  


# ----------------------------------------------------------------------
# >>> Initialize empty complex 2D-array of given                  ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/05/08  - adapted from L.Laguarda
#
# Desc
#  - row-majored : order='C'
# ----------------------------------------------------------------------

def init_2Dcmx_empty( rows, cols, precision ):

   if precision == 8:
      return np.empty( (rows,cols), dtype='c8', order='C' )

   else:
      return np.empty( (rows,cols), dtype='c16', order='C' )