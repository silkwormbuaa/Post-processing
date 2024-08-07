# -*- coding: utf-8 -*-
'''
@File    :   io_vtk.py
@Time    :   2024/08/07 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import vtk
import numpy                         as     np
from   vtkmodules.util               import numpy_support
from   vtkmodules.vtkCommonDataModel import vtkRectilinearGrid 


# ----------------------------------------------------------------------
# >>> create_3d_vtkRectilinearGrid                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def create_3d_vtkRectilinearGrid(x, y, z):
    
    """
    create a 3d vtkRectilinearGrid object with x, y, z
    
    return: vtkRectilinearGrid object with points defined.
    """
    
    grid = vtkRectilinearGrid()
    grid.SetDimensions( len(x), len(y), len(z) )
    
    grid.SetXCoordinates( numpy_support.numpy_to_vtk(x) )
    grid.SetYCoordinates( numpy_support.numpy_to_vtk(y) )
    grid.SetZCoordinates( numpy_support.numpy_to_vtk(z) )
    
    return grid

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def add_var_vtkRectilinearGrid(grid:vtkRectilinearGrid, var_name, var_value):
    
    """
    add a variable to the vtkRectilinearGrid object
    
    grid: vtkRectilinearGrid object
    var_name: string of the variable name
    var_value: numpy array with column major order (z,y,x)
    """
    
    var = numpy_support.numpy_to_vtk( np.ravel(var_value) )
    
    var.SetName(var_name)
    
    var.SetNumberOfComponents(1)
    var.SetNumberOfTuples( grid.GetNumberOfCells() )
    
    grid.GetCellData().AddArray(var)
    
    return grid

# ----------------------------------------------------------------------
# >>> create_multiblock_dataset                                 (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def create_multiblock_dataset(blocks):
    
    """
    Create a vtkMultiBlockDataSet from a list of structured grid blocks.

    blocks: List of vtkStructuredGrid objects.
    return: vtkMultiBlockDataSet object.
    """
    
    multi_block_dataset = vtk.vtkMultiBlockDataSet()
    
    for i, block in enumerate(blocks):
        multi_block_dataset.SetBlock(i, block)
    
    return multi_block_dataset

# ----------------------------------------------------------------------
# >>> write_vtm_file                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def write_vtm_file(filename, multiblock_dataset):
    
    """
    write the multiblock dataset to a vtm file.
    
    filename: name of the file
    multiblock_dataset: vtkMultiBlockDataSet object
    """

    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(multiblock_dataset)
    writer.Write()
    
    
# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():
    
    os.chdir('/home/wencanwu/test/vtk/')

    x = np.linspace(0, 1, 10)
    y = np.linspace(0, 2, 20)
    z = np.linspace(0, 3, 30)
    
    grid = create_3d_vtkRectilinearGrid(x, y, z)
    
    u = np.zeros( (len(x)-1, len(y)-1, len(z)-1) )
    
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            for k in range(len(z)-1):
                u[i,j,k] = x[i]+y[j]+z[k]
        
    grid = add_var_vtkRectilinearGrid(grid, 'u', u.T)
    
    print(grid.GetNumberOfCells())
    print(grid.GetNumberOfPoints())
    print(grid.GetCellData().GetArray('u').GetNumberOfTuples())
    print(grid.GetCellData().GetArray('u').GetNumberOfComponents())
    print(grid.GetCellData().GetArray('u').GetName())
    
    dataset = create_multiblock_dataset([grid])
    write_vtm_file('test.vtm', dataset)
    

# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/08/07  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
