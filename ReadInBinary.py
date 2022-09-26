# -*- coding: utf-8 -*-
'''
@File    :   ReadInBinary.py
@Time    :   2022/09/14 12:19:49
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import numpy as np

# ------------------------------------------------------------------------
# READ SIGNED INTEGER(S) FROM BYTE STREAM                           (  1 )
# L. Laguarda
#
# March, 2021  - created
# ------------------------------------------------------------------------
def read_int_bin( byte_stream, precision ):

   # Get number of integers in byte_stream
   n = len( byte_stream ) % precision

   if len( byte_stream ) % precision == 0:
      n = int(len( byte_stream )/precision)
   else:
      fatal_error_mp( 'Inconsistent byte stream' )

   # Read data (if output is an array, it is in column-major order)
   if   precision == 2:

      if n > 1:
         return np.array( struct.unpack_from( n*'h',
                          byte_stream, offset=0), dtype='i2', order='F' )
      else:
         return struct.unpack_from( 'h', byte_stream, offset=0)[0]

   elif precision == 4:

      if n > 1:
         return np.array( struct.unpack_from( n*'i',
                          byte_stream, offset=0), dtype='i4', order='F' )
      else:
         return struct.unpack_from( 'i', byte_stream, offset=0)[0]

   elif precision == 8:

      if n > 1:
         return np.array( struct.unpack_from( n*'d',
                          byte_stream, offset=0), dtype='i8', order='F' )
      else:
         return struct.unpack_from( 'q', byte_stream, offset=0)[0]

   else:
      fatal_error_mp( 'Float precision not supported' )

# ------------------------------------------------------------------------
# READ REAL NUMBER(S) FROM BYTE STREAM                              (  2 )
# L. Laguarda
#
# March, 2021  - created
# ------------------------------------------------------------------------
def read_flt_bin( byte_stream, precision ):

   # Get number of floats in byte_stream
   n = len( byte_stream ) % precision

   if len( byte_stream ) % precision == 0:
      n = int(len( byte_stream )/precision)
   else:
      fatal_error_mp( 'Inconsistent byte stream' )

   # Read data (if output is an array, it is in column-major order)
   if   precision == 4:

      if n > 1:
         return np.array( struct.unpack_from( n*'f',
                          byte_stream, offset=0), dtype='f4', order='F' )
      else:
         return struct.unpack_from( 'f', byte_stream, offset=0)[0]

   elif precision == 8:

      if n > 1:
         return np.array( struct.unpack_from( n*'d',
                          byte_stream, offset=0), dtype='f8', order='F' )
      else:
         return struct.unpack_from( 'd', byte_stream, offset=0)[0]


   else:
      fatal_error_mp( 'Float precision not supported' )

# ------------------------------------------------------------------------
# READ BLOCK SNAPSHOT HEADER - BINARY                               (  2 )
# L. Laguarda
#
# September, 2021  - created
# ------------------------------------------------------------------------
   def read_header_binary( self, file ):

      # Book-keeping
      pos = 0

      # Default float size
      sin = 4
      slg = 4
      sfl = 8

      # Header format
      hformat = read_int_bin( file.read(sin), sin )

      # Floating point precision format
      kind    = read_int_bin( file.read(sin), sin )

      # Time step
      itstep  = read_int_bin( file.read(sin), sin )

      # Number of species
      nspec   = read_int_bin( file.read(sin), sin )
      pos += 4*sin

      # Simulation time
      itime   = read_flt_bin( file.read(sfl), sfl )
      pos += sfl

      # Content
      buf_log = read_log_bin( file.read(8*slg), slg )
      pos += 8*slg

      # Assign to class variables
      self.kind          = kind
      self.itstep        = itstep
      self.itime         = itime
      self.n_species     = nspec
      self.snap_lean     = buf_log[0]
      self.compressible  = buf_log[1]
      self.snap_with_gx  = buf_log[2]
      self.snap_with_tp  = buf_log[3]
      self.snap_with_vp  = buf_log[4]
      self.snap_with_cp  = buf_log[5]
      self.snap_with_mu  = buf_log[6]
      self.snap_with_wd  = buf_log[7]
