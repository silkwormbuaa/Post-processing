# -*- coding: utf-8 -*-
'''
@File    :   io_directories.py
@Time    :   2024/03/04 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   Module of input/output directories names
'''


import os
from   .tools            import get_filelist

class Directories:
    
    def __init__(self, case_dir):
        
        """
        case_dir: str, path to the case directory
        
        generate static directories under current case directory
        """
        
        # -- first level directories
        self.case_dir = os.path.abspath( case_dir )
        self.res_dir  = os.path.join( self.case_dir, 'results' )
        self.set_dir  = os.path.join( self.case_dir, 'setup' )
        self.tec_dir  = os.path.join( self.case_dir, 'tec' )
        self.snp_dir  = os.path.join( self.case_dir, 'snapshots' )
        self.prb_dir  = os.path.join( self.case_dir, 'probes' )
        self.sup_dir  = os.path.join( self.case_dir, 'supplements' )
        self.pos_dir  = os.path.join( self.case_dir, 'postprocess' )

        # -- second level directories
        
        # results directories items
        
        self.statistics = os.path.join( self.res_dir, 'statistics.bin' )
        self.solution   = os.path.join( self.res_dir, 'results.bin' )
        self.grid       = os.path.join( self.res_dir, 'inca_grid.bin' )
          
        
        # tecplot szplt files 
        
        self.tec_stat = os.path.join( self.tec_dir, 'TP_stat' )
        
        # snapshots
        
        # snapshots directories
        
        # probe files
        
        # supplement files
        
        self.cc_setup = os.path.join( self.sup_dir, 'cutcells_setup.dat' )
        self.cc_BSP_setup = os.path.join( self.sup_dir, 'cutcells_BSP_setup.dat' )
        self.wall_dist = os.path.join( self.sup_dir, 'wall_dist' )
        self.case_para_file = os.path.join( self.sup_dir, 'case_parameters' )

        # post processing directoris
        
        self.pp_statistics   = os.path.join( self.pos_dir, 'statistics' )
        self.pp_probes       = os.path.join( self.pos_dir, 'probes' )
        self.pp_solution     = os.path.join( self.pos_dir, 'solution' )
        self.pp_snapshots    = os.path.join( self.pos_dir, 'snapshots' )
        self.pp_DMD          = os.path.join( self.pos_dir, 'DMD' )
        
        self.pp_wall_proj    = os.path.join( self.pp_statistics, 'wall_projection' )
        
# ----------------------------------------------------------------------
# >>> Dynamic files                                            (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/04  - created
#
# Desc
#
# ----------------------------------------------------------------------

    # results
    @property
    def stl( self ):
        return get_filelist( self.res_dir, '.stl' )[0]

    # snapshots
    @property
    def snaps( self ):
        items = os.listdir( self.snp_dir )
        folders = [ item for item in items if os.path.isdir( os.path.join( self.snp_dir, item ) ) ]
        return os.listdir( folders )
    
        
# ----------------------------------------------------------------------
# >>> Testing section                                           ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/04  - created
#
# Desc
#
# ----------------------------------------------------------------------

def Testing():

    case = Directories('.')
    print(case.res_dir)


# ----------------------------------------------------------------------
# >>> Main: for test and debugging                              ( -1 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/04  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":

    Testing()
