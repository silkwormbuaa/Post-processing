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
from   pathlib           import Path
from   .tools            import get_filelist

class Directories:
    
    def __init__(self, case_dir):
        
        """
        case_dir: str, path to the case directory\n
        generate static directories under current case directory
        """
        
    # --- first level directories
        self.case_dir = os.path.abspath( case_dir )
        self.res_dir  = os.path.join( self.case_dir, 'results' )
        self.set_dir  = os.path.join( self.case_dir, 'setup' )
        self.tec_dir  = os.path.join( self.case_dir, 'tec' )
        self.snp_dir  = os.path.join( self.case_dir, 'snapshots' )
        self.prb_dir  = os.path.join( self.case_dir, 'probes' )
        self.sup_dir  = os.path.join( self.case_dir, 'supplements' )
        self.pos_dir  = os.path.join( self.case_dir, 'postprocess' )

        # -- second level directories
        
        # results directory items
        
        self.statistics = os.path.join( self.res_dir, 'statistics.bin' )
        self.solution   = os.path.join( self.res_dir, 'results.bin' )
        self.grid       = os.path.join( self.res_dir, 'inca_grid.bin' )
        
        # setup directory items  
        
        self.set_prb = os.path.join( self.set_dir, 'inca_probes.inp' )
        
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

# - post processing directoris
        
        self.pp_statistics   = os.path.join( self.pos_dir, 'statistics' )
        self.pp_probes       = os.path.join( self.pos_dir, 'probes' )
        self.pp_solution     = os.path.join( self.pos_dir, 'solution' )
        self.pp_snapshots    = os.path.join( self.pos_dir, 'snapshots' )
        self.pp_DMD          = os.path.join( self.pos_dir, 'DMD' )
        
# --- pp_statistics subdirectories
        
        self.pp_wall_proj     = os.path.join( self.pp_statistics, 'wall_projection' )
        self.pp_profile_array = os.path.join( self.pp_statistics, 'profile_array' )
        self.pp_integral      = os.path.join( self.pp_statistics, 'integral' )
        
# --- pp_probes subdirectories
        
        self.pp_psd_ridge  = os.path.join( self.pp_probes, 'psd_ridge'  )
        self.pp_psd_valley = os.path.join( self.pp_probes, 'psd_valley' )
        self.pp_psd_others = os.path.join( self.pp_probes, 'psd_others' )
        self.pp_pre_ridge  = os.path.join( self.pp_probes, 'pre_ridge'  )
        self.pp_pre_valley = os.path.join( self.pp_probes, 'pre_valley' )
      
# --- pp_snapshots subdirectories

        self.pp_snp_fricpv  = os.path.join( self.pp_snapshots, 'fricpv'  )
        self.pp_snp_fricprj = os.path.join( self.pp_snapshots, 'fricprj' )
        self.pp_snp_pfmax   = os.path.join( self.pp_snapshots, 'pfmax'   )
        
        self.pp_signals    = os.path.join( self.pp_probes, 'signals'    )

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
    
    # probes list
    @property
    def probes( self ):
        return get_filelist( self.prb_dir, '.dat' )
    

# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/21  - created
#
# Desc
#
# ----------------------------------------------------------------------

def create_folder( path ):
    
    try:
        Path(path).mkdir(parents=True, exist_ok=True)
        print(f"Folder created at path: {path}")
    except FileExistsError:
        print(f"Folder already exists at path: {path}")
    
    return path

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
