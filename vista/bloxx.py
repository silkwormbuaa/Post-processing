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
import pandas            as     pd

from   .tools            import get_filelist

class Grid_bloxx:
    def __init__(self, file_path):
        self.variables = self._read_grid_file(file_path)
        self.filename  = file_path.split('/')[-1]

    def _read_grid_file(self, file_path):
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


# ----------------------------------------------------------------------
# >>> Write Grid_bloxx file                                      (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/11/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def write_to_file(self, output_file_path=None):
        
        """
        output_file_path: output file path
        Write the grid to a file.
        """

        output_file_path = output_file_path or self.filename
        
        with open(output_file_path, 'w') as file:
            file.write("! Automatically generated INCA grid file.\n\n")
            file.write("&GRID\n")

            for variable, values in self.variables.items():
                # Write variable name
                file.write(f" {variable}={values}\n")

            file.write(" /\n")


# ----------------------------------------------------------------------
# >>> Define attributes                                           (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2023/11/23  - created
#
# Desc
#
# ----------------------------------------------------------------------

    @property
    def LX(self):
        LX_str = self.variables['LX']
        numbers = LX_str.split(',')
        return tuple( float(number) for number in numbers if number.strip() )

    @property
    def LY(self):
        LY_str = self.variables['LY']
        numbers = LY_str.split(',')
        return tuple( float(number) for number in numbers if number.strip() )
    
    @property
    def LZ(self):
        LZ_str = self.variables['LZ']
        numbers = LZ_str.split(',')
        return tuple( float(number) for number in numbers if number.strip() )
    
    @property
    def NX(self):
        NX_str = self.variables['NX']
        return int( NX_str.strip(',').strip() )
    
    @property
    def NY(self):
        NY_str = self.variables['NY']
        return int( NY_str.strip(',').strip() )
    
    @property
    def NZ(self):
        NZ_str = self.variables['NZ']
        return int( NZ_str.strip(',').strip() )
    
    @property
    def BX(self):
        BX1_str = self.variables['BX1'].strip(',').strip().strip('"')
        BX2_str = self.variables['BX2'].strip(',').strip().strip('"')
        return (BX1_str, BX2_str)
    
    @property
    def BY(self):
        BY1_str = self.variables['BY1'].strip(',').strip().strip('"')
        BY2_str = self.variables['BY2'].strip(',').strip().strip('"')
        return (BY1_str, BY2_str)
    
    @property
    def BZ(self):
        BZ1_str = self.variables['BZ1'].strip(',').strip().strip('"')
        BZ2_str = self.variables['BZ2'].strip(',').strip().strip('"')
        return (BZ1_str, BZ2_str)

    @property
    def BC(self):
        return self.BX + self.BY + self.BZ

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/30  - created
#
# Desc
#
# ----------------------------------------------------------------------

class Mesh_bloxx:
    
    def __init__(self, folder):
        
        self.folder =  folder
        self.grid_files = get_filelist(folder, 'inca_grid')
        
        self.grids = []
        for file in self.grid_files:
            grid = Grid_bloxx(file)
            self.grids.append(grid)
        
        
# ----------------------------------------------------------------------
# >>> check_boundary_conditions                                  (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/30  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def check_boundary_conditions( self ):
        
        """
        check how many times which boundary conditions are used
        """
        
        # Initialize a dictionary to count boundary conditions
        boundary_conditions = {
            'CYC': 0,
            'DF_INFLOW': 0,
            'RI_INFLOW': 0,
            'OUTFLOW': 0,
            'DUMMY':0
        }
        
        for grid in self.grids:
            
            # count boundary conditions
            
            for condition in boundary_conditions.keys():
                
                boundary_conditions[condition] += grid.BC.count(condition)
            
        print(boundary_conditions)
    

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/01/30  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def sort_grids(self, output_folder):
        
        """
        output_folder: output folder path
        resort grids in the output folder
        """
        
        # create output folder if not exist
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        
        lx1 = [grid.LX[1] for grid in self.grids]
        ly1 = [grid.LY[1] for grid in self.grids]
        lz1 = [grid.LZ[1] for grid in self.grids]
        old_filename = [grid.filename for grid in self.grids]
        
        df = pd.DataFrame({'lx1': lx1, 'ly1': ly1, 'lz1': lz1, 
                           'old_filename': old_filename})
        
        # sort grids and reset index
        df = df.sort_values(by=['lx1', 'ly1', 'lz1'])
        df = df.reset_index(drop=True)
        
        print(df)
        
        # generate a dict between filename and index
        new_index = {filename: index for filename, index in zip(df['old_filename'], df.index+1)}
        
        # change grid filename and write to file
        for grid in self.grids:
            
            grid.filename = f"inca_grid_{new_index[grid.filename]:06d}.inp"
            grid.write_to_file(output_folder + '/' + grid.filename)
        
        print(f"sorted grids are output to {output_folder}")
            
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
    file_path = '/home/wencanwu/my_simulation/STBLI_mid_Re/grid_experiment/sorted_grids/'
    
    os.chdir('/home/wencanwu/my_simulation/STBLI_mid_Re/grid_experiment/grid')
    
    mesh = Mesh_bloxx(file_path)
    
    mesh.sort_grids('../sorted_grids/')


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
