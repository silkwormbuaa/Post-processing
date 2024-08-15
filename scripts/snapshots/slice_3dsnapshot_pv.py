#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   slice_3dsnapshot_pv.py
@Time    :   2024/08/15 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Get the snapshot using pyvista
'''

import os
import sys
import pyvista           as     pv

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.timer       import timer
from   vista.log         import Logger
sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

snapshotfile = 'snapshot.bin'
gridfile = 'grid.bin'

slic_vector = [0,1,0]
slic_origin = [0,0,0]

vars = ['u','v','w','p','T']

# - read grid file

G = GridData( gridfile )
G.read_grid()

snapshot = Snapshot( snapshotfile )
snapshot.grid3d = G
snapshot.read_snapshot()

# - slicing vector

dataset = snapshot.create_vtk_multiblock( vars )

dataset.set_active_scalars('u')

slice = dataset.slice(normal=slic_vector, origin=slic_origin)

p = pv.Plotter()

p.add_mesh(slice, opacity=1.0, show_scalar_bar=True)

p.show()