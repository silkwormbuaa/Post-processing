#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   read_binary.py
@Time    :   2023/02/01 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Basic functions of reading binary data.
'''

# from multiprocessing.sharedctypes import Value

import struct

import numpy             as np

# ----------------------------------------------------------------------
# >>> read_int_bin                                               ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu   : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/01  - created
#
# 2021/03     - originally from L. Laguarda  
#
# Desc
#
# - read signed integer(s) from byte stream
#
# ----------------------------------------------------------------------


def read_int_bin( byte_stream, precision ):


    # Get number of integers in byte_stream

    n = len( byte_stream ) % precision

    if len( byte_stream ) % precision == 0:

        n = int(len( byte_stream )/precision)

    else:
        raise ValueError('Inconsistent byte stream')

    # Read data (if output is an array, it is in column-major order)

    if precision == 2:


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
        raise ValueError('Int precision not supported')

# ----------------------------------------------------------------------
# >>> read_flt_bin                                                ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/01  - created
#
# Desc
#
# - read real number(s) from byte stream
#
# ----------------------------------------------------------------------


def read_flt_bin( byte_stream, precision ):


    # Get number of floats in byte_stream

    n = len( byte_stream ) % precision


    if len( byte_stream ) % precision == 0:

        n = int(len( byte_stream )/precision)

    else:
        raise ValueError('Inconsistent byte stream')


    # Read data (if output is an array, it is in column-major order)

    if precision == 4:


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
        raise ValueError('Float precision not supported')

# ----------------------------------------------------------------------
# >>> read_log_bin                                                ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/01  - created
#
# Desc
#
# - read logical(s) from byte stream
#
# ----------------------------------------------------------------------


def read_log_bin( byte_stream, precision ):


    # Get number of floats in byte_stream

    n = len( byte_stream ) % precision


    if len( byte_stream ) % precision == 0:

        n = int(len( byte_stream )/precision)

    else:
        raise ValueError('Inconsistent byte stream')


    # Read data (if output is an array, it is in column-major order)

    if precision == 1 or precision == 4:

        if n > 1:

#            buf_log = init_1Dlog_empty( n )
            buf_log = [None]*n

            ib = 0

            fb = precision


            for i in range(0,n):

                bs = byte_stream[ib:fb]

                buf_log[i] = struct.unpack_from( '?', bs, offset=0 )[0]

                ib += precision

                fb += precision

            return buf_log


        else:

            return struct.unpack_from( '?', byte_stream, offset=0)[0]


    else:
        raise ValueError('Logic precision not supported')

# ----------------------------------------------------------------------
# >>> read_chr_bin                                                ( 4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/02/01  - created
#
# Desc
#
# - read character(s) from byte stream
#
# ----------------------------------------------------------------------


def read_chr_bin( byte_stream ):

    # Get number of floats in byte_stream

    n = len( byte_stream ) 

    # Read data (if output is an array, it is in column-major order)
    return struct.unpack_from( '%ds'%n, byte_stream
                             , offset = 0 )[0].decode('utf-8')


# ----------------------------------------------------------------------
# >>> read_bin_3Dflt                                            ( 4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/04/24  - created
#
# Desc
#
# ----------------------------------------------------------------------

def read_bin_3Dflt( pos, file, N1, N2, N3, kind ):
    
    # Get number of floats that should be read
    
    n = N1 * N2 * N3
    
    # move the file pointer to where reading starts
    
    file.seek( pos )
    
    flt_buf = read_flt_bin( file.read( n*kind ), kind )
    
    return flt_buf