#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   convert_fig.py
@Time    :   2023/10/17 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :       
    Now can only convert from bitmap to vector image, cannot the
    other way around.
    
    ! Suggest generate all format figures at the beginning. When there is really 
      need to convert figures, use linux command 'convert'
'''


import os

import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.tools       import get_filelist
from   vista.tools       import convert_image_format


# receive parameters from user

if len( sys.argv ) not in [3,5]:
    print(">> python3 convert_fig.py <input_format> <output_format> [w] [h]")
    sys.exit(1)

input_format  = sys.argv[1]
output_format = sys.argv[2].lower()

output_size   = None

if len( sys.argv ) == 5:
    try:
        width  = int(sys.argv[3])
        height = int(sys.argv[4]) 
        output_size = ( width, height )
    except ValueError:
        print("width and height must be integers!")
        sys.exit(1)


# enter the output folder

imgpath = os.getcwd()

output_dir = imgpath + '/'+ 'converted_' + output_format

if not os.path.exists(output_dir): 
    os.mkdir( output_dir )
    print(f"Created directory {output_dir}.\n")

os.chdir(output_dir)

# find all images with input format

image_files = get_filelist( imgpath, '.'+input_format )

for image_file in image_files:
    
    input_dir, input_name_with_ext = os.path.split( image_file) 
    input_name, _ = os.path.splitext( input_name_with_ext )
    output_file = input_name + '.' + output_format
    
    convert_image_format( image_file, output_file, output_format, 
                          outputsize= output_size )
