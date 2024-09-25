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

snapshotfile = '/home/wencanwu/test/220927/snapshots/snapshot_01327770/snapshot_Z_001.bin'
gridfile = '/home/wencanwu/test/220927/results/inca_grid.bin'

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

# pass snapshot to pyvista

vars = ['u','v','w','p','T','grad_rho','grad_p']
rescale = [0.0,.0,.0,1.0,1.0,1.0]
# pv.global_theme.font.color = 'white'
# pv.global_theme.font.size  = 10
dataset = pv.MultiBlock( snap.create_vtk_multiblock( vars, mode='symmetry', rescale=rescale ) )

dataset.set_active_scalars('grad_rho')

p = pv.Plotter()
print(dataset)
print(dataset.bounds)
p.add_mesh(dataset, opacity=1.0, clim=[0.0,0.6],above_color='red', below_color='blue', show_scalar_bar=True)
# p.set_background('black')

p.show_bounds(  dataset, 
              bounds=[-116,120,0,80,0,0],
              axes_ranges=[-116,120,0,80,0,0],
               ticks='both', show_zaxis=False,
               xtitle=r'$(x-x_{imp})/\delta_0$',
               ytitle=r'$y/\delta_0$',
               font_size=10)

p.view_vector([0.0,0.0,1.0],viewup=[0.0,1.0,0.0])

p.show_axes()
p.show()



