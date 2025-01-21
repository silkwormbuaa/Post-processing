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
import shutil
import numpy             as     np
import pandas            as     pd
from   copy              import deepcopy

from   .tools            import get_filelist
from   .directories      import create_folder

class Grid_bloxx:

    """
    file_path: path to the grid file, if None, create an empty grid \n
    """

    def __init__(self, file_path=None):
        
        if file_path is None:
            self.filename  = ''
            self.variables = {}

            variables_to_set = {
                'LX       ': '0.0, 1.0 ,',
                'LY       ': '0.0, 1.0 ,',
                'LZ       ': '0.0, 1.0 ,',
                'NX       ': '1,',
                'NY       ': '1,',
                'NZ       ': '1,',
                'TIMES_PI ': ' .F. ,',
                'SHAPEX   ': '"HOMO" ,',
                'PARAMX   ': '0.00000000e+00  , 0.00000000e+00  , 0.00000000e+00  ,',
                'SHAPEY   ': '"HOMO" ,',
                'PARAMY   ': '0.00000000e+00  , 0.00000000e+00  , 0.00000000e+00  ,',
                'SHAPEZ   ': '"HOMO" ,',
                'PARAMZ   ': '0.00000000e+00  , 0.00000000e+00  , 0.00000000e+00  ,',
                'BX1      ': '"DUMMY" ,',
                'FLUIDX1  ': '"DUMMY" ,',
                'VALX1    ': '0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,',
                'PX1      ': '0               , 0               , 0               ,',
                'BX2      ': '"DUMMY" ,',
                'FLUIDX2  ': '"DUMMY" ,',
                'VALX2    ': '0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,',
                'PX2      ': '0               , 0               , 0               ,',
                'BY1      ': '"DUMMY" ,',
                'FLUIDY1  ': '"DUMMY" ,',
                'VALY1    ': '0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,',
                'PY1      ': '0               , 0               , 0               ,',
                'BY2      ': '"DUMMY" ,',
                'FLUIDY2  ': '"DUMMY" ,',
                'VALY2    ': '0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,',
                'PY2      ': '0               , 0               , 0               ,',
                'BZ1      ': '"DUMMY" ,',
                'FLUIDZ1  ': '"DUMMY" ,',
                'VALZ1    ': '0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,',
                'PZ1      ': '0               , 0               , 0               ,',
                'BZ2      ': '"DUMMY" ,',
                'FLUIDZ2  ': '"DUMMY" ,',
                'VALZ2    ': '0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,',
                'PZ2      ': '0               , 0               , 0               ,',
            }
            
            for key, value in variables_to_set.items():
                self.variables[key.strip()] = value
        
        else:
            self.filename  = file_path.split('/')[-1]
            self.variables = self._read_grid_file(file_path)

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
        output_file_path: output file path \n
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
# >>> class of Mesh_bloxx                                          ( 0 )
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
    
    def __init__(self, folder=None):
        
        """
        folder: folder path containing grid files
        
        read each grid files as Grid_bloxx object and store them in self.grids
        """
        
        if folder is None:
            
            self.folder     = ''
            self.grid_files = []
            self.grids      = []
    
        else:
        
            self.folder     = folder
            self.grid_files = get_filelist( folder, 'inca_grid' )
            
            self.grids = []
            
            for file in self.grid_files:
                
                grid = Grid_bloxx(file)
                self.grids.append(grid)
        
        
# ----------------------------------------------------------------------
# >>> check_boundary_conditions                                  ( 1 )
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
# >>> sort_grids                                                ( 2 )
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

    def sort_grids(self):
        
        """
        resort grids in self.grids based on (lx1,ly1,lz1)
        """
        
        lx1 = [grid.LX[1] for grid in self.grids]
        ly1 = [grid.LY[1] for grid in self.grids]
        lz1 = [grid.LZ[1] for grid in self.grids]
        old_filename = [grid.filename for grid in self.grids]
        
        df = pd.DataFrame({'lx1': lx1, 'ly1': ly1, 'lz1': lz1, 'grids': self.grids,
                           'old_filename': old_filename})
        
        # sort grids and reset index
        
        df = df.sort_values(by=['lx1', 'ly1', 'lz1'])
        df = df.reset_index(drop=True)
        
        for i, grid in enumerate(df['grids']):
            grid.filename = f"inca_grid_{i+1:06d}.inp"
        
        
        # add new filenames to the dataframe
        
        new_filename = [grid.filename for grid in df['grids']]
        df['new_filename'] = new_filename
        
        print( df )
        
        print(f"Grids are sorted based on (x,y,z) !")


# ----------------------------------------------------------------------
# >>> save grids                                                   ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/02/14  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def save_grid(self, output_folder):
        
        """
        output_folder: output folder path
        """
        
        if len(self.grids) == 0:
            
            print("\033[93m No grid to save \033[0m")
            return
        
        else:

            # create output folder if not exist
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            
            for grid in self.grids:
                grid.write_to_file(output_folder.rstrip('/') + '/' + grid.filename)

            print(f"\033[92mGrids are saved to {output_folder} \033[0m")
            
# ----------------------------------------------------------------------
# >>> Select blocks in a given box                                ( 4 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/02/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

    def select_blocks(self, box):
        
        """
        box: a tuple of 6 floats (x1, y1, z1, x2, y2, z2)\n
        select blocks in the box
        """
        
        selected_grids = []
        
        for grid in self.grids:
            
            if (grid.LX[0] >= box[0] and grid.LX[1] <= box[3] and
                grid.LY[0] >= box[1] and grid.LY[1] <= box[4] and
                grid.LZ[0] >= box[2] and grid.LZ[1] <= box[5]):
                
                selected_grids.append(grid)
        
        return selected_grids
    
    

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

    pass

# Example usages:

# =============================================================================

# -- sort grids based on (lx1,ly1,lz1)

#     file_path = '/home/wencanwu/test/bloxx_test/grid'
    
#     os.chdir(file_path.split('grid')[0])
    
#     mesh = Mesh_bloxx(file_path)
    
# #    mesh.grids = mesh.select_blocks((-120.0,-0.30,-11.0,120.0,100.0,11.0))
    
#     mesh.sort_grids()

#     mesh.save_grid( create_folder('./new_grid') )
    
#     grid = mesh.grids[0]
    
#     for variable in grid.variables:
#         print(variable, grid.variables[variable])


# =============================================================================

# # -- add two layers of new blocks at the bottom of the domain

#     new_grids = []
#     i_new     = len( mesh.grids )
    
#     for grid in mesh.grids:
        
#         if grid.LY[0] == -0.7524:
            
#             i_new += 1
#             grid_new = deepcopy(grid)
#             grid_new.variables['LY'] = ' -1.0032, -0.7524'
#             grid_new.filename = f'inca_grid_{i_new:06d}.inp'
#             new_grids.append(grid_new)
            
#             i_new += 1
#             grid_new = deepcopy(grid)
#             grid_new.variables['LY'] = ' -1.2540, -1.0032'
#             grid_new.filename = f'inca_grid_{i_new:06d}.inp'
#             new_grids.append(grid_new)
        
#     mesh.grids += new_grids

# =============================================================================

# # -- modifiy boundary conditions

#     for grid in mesh.grids:
        
#         if grid.LY[0] == 0.0:
#             grid.variables['BY1'] = 'AWALL     '
#             print(grid.BC)

# =============================================================================



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
