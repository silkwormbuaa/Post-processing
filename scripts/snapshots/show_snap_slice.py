    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   shock_tracking.py
@Time    :   2024/08/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   read in a snapshot and show it with pyvista
'''


import os
import sys
import pyvista           as     pv

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.snapshot    import Snapshot
from   vista.grid        import GridData
from   vista.timer       import timer
#from   vista.log         import Logger
#sys.stdout = Logger( os.path.basename(__file__) )

# =============================================================================

snapshotfile = '/media/wencanwu/Seagate Expansion Drive1/temp/231124/snapshots/snapshot_02920790/snapshot_Y_003.bin'
gridfile = '/media/wencanwu/Seagate Expansion Drive1/temp/231124/results/inca_grid.bin'

# =============================================================================

# - read grid file

grd = GridData( gridfile )
grd.read_grid()

# - read snapshot file

snap = Snapshot( snapshotfile )
snap.grid3d = grd
snap.read_snapshot( )

# compute gradients in the snapshot

snap.compute_gradients( grads=['grad_rho','grad_p'] )

print( max([max(bl.df['grad_p']) for bl in snap.snap_data]) )

# pass snapshot to pyvista

vars = ['u','v','w','p','T','grad_rho','grad_p']
dataset = pv.MultiBlock(snap.create_vtk_multiblock( vars ))

print( type(dataset) )

dataset.set_active_scalars('grad_rho')

# grad_p [0,50000]

p = pv.Plotter()
print(dataset)
p.add_mesh(dataset, opacity=1.0, clim=[0,0.5], above_color='red', below_color='blue', show_scalar_bar=True)
p.show_axes()
axes_actor = pv.CubeAxesActor(p.camera)
axes_actor.bounds = dataset.bounds
actor, property = p.add_actor( axes_actor)
p.show()



