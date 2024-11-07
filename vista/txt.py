# -*- coding: utf-8 -*-
'''
@File    :   txt.py
@Time    :   2022/11/29 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Methods used for manipulating texts.
'''

import os

# ----------------------------------------------------------------------
# >>> Replace text                                              ( 1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/11/29  - created
#
# Desc
# 
# - replace the text 
#
# ----------------------------------------------------------------------

def txt_replace( file, in_str, out_str ):

    with open( file, 'r' ) as f:
        
        lines = f.readlines()
        
        content = ''
        
        for line in lines:
            
            line = line.replace( in_str, out_str )
            
            content = content + line
        
        print( content )

    with open( file, 'w' ) as f:
        
        f.write( content )
        

# ----------------------------------------------------------------------
# >>> tail a text file                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/11/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def txt_tail( filename, lines=10 ):

    """
    filename: ascii file name
    lines: number of lines to read from the end
    
    return: the last few lines of the file
    """

    with open(filename, "rb") as f:
        
        f.seek(0, os.SEEK_END)
        buffer = bytearray()
        pointer = f.tell()
        
        # - read the byte stream from the end, stored in an inversed buffer
        
        while pointer >= 0 and lines > 0:
            f.seek(pointer)
            pointer -= 1
            new_byte = f.read(1)
            if new_byte == b'\n':
                lines -= 1
            buffer.extend(new_byte)
        
        # - flip the buffer and decode 
        
        return buffer[::-1].decode()


# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/11/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():
    
    
    filename = '/media/wencanwu/Seagate Expansion Drive1/temp/test/run.out'
    
    message = txt_tail( filename, 10 )
    
    keyword = "WARNING"
    
    if keyword in message:
        
        print( "Keyword found!" )



# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/11/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
