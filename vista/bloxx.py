# -*- coding: utf-8 -*-
'''
@File    :   bloxx.py
@Time    :   2023/11/22 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os

from   .tools            import get_filelist

class Grid_bloxx:
    def __init__(self, file_path):
        self.variables = self._parse_inca_file(file_path)

    def _parse_inca_file(self, file_path):
        variables = {}

        with open(file_path, 'r') as file:
            lines = file.readlines()

            for line in lines:
                # split from the first '='
                
                if line.strip() and not (line.strip().startswith(("!","&","/"))):

                    key, value = line.strip().split("=",1) # split once
                    
                    # remove extra space and comma
                    
                    key = key.strip()
                    
                    variables[key] = value

        return variables
    

    def write_to_file(self, output_file_path):
        with open(output_file_path, 'w') as file:
            file.write("! Automatically generated INCA grid file.\n\n")
            file.write("&GRID\n")

            for variable, values in self.variables.items():
                # Write variable name
                file.write(f" {variable}={values}\n")

            file.write(" /\n")




# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/11/22  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    # Example usage:
    file_path = '/home/wencanwu/my_simulation/STBLI_mid_Re/grid_ascii/'
    files = get_filelist(file_path)
    
    os.chdir(file_path)
    for i, file in enumerate(files):
        
        inca_instance = Grid_bloxx(file)
        print(i,len(inca_instance.variables))


    # Modify the instance as needed
    # ...

    # Write the modified instance to a new file
#    output_file_path = '/home/wencanwu/my_simulation/STBLI_mid_Re/grid_test.inp'
#    inca_instance.write_to_file(output_file_path)

# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/11/22  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
