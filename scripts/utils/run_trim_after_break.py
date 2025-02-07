#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   vista_trim_after_break.py
@Time    :   2022/12/07 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   - Cleaning a INCA case after accidental break
             - Trim 
              a. probes
              b. snapshots directories and snapshot_index
              c. force on each block
'''

import os
import sys

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.tools       import get_filelist

# ----------------------------------------------------------------------
# >>> Control Panel                                               ( 0 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/12/07  - created
#
# Desc
#
# - trim_ibforce need backup time from trim_snapshots
# ----------------------------------------------------------------------

# input the step of your backup (where you want to continue computation)

backupstep = int( 1340310 )

trim_snapshots = True
trim_probe     = True
trim_ibforce   = True


folder = os.getcwd()

#folder = '/home/wencanwu/testdir/case'

print( folder )

# --- Fail-safe measures

if folder[0:35] == '/home/wencanwu/Code/Post-processing':
    
    raise KeyError("The working path is wrong!")

else: 

    print( 'start cleanning ...\n')

# ----------------------------------------------------------------------
# >>> Cleaning Snapshots                                          ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/12/07  - created
#
# Desc
#
# - trim the snapshot_index.dat and also get the backup time
# - delete the snapshots directories after backup step 
#
# ----------------------------------------------------------------------

if trim_snapshots:
    
    os.chdir( folder )

    os.chdir( 'snapshots' )
    
    print( os.getcwd() )
    
# --- trim snapshot_index.dat 

    with open( 'snapshot_index.dat', 'r+b' ) as f:
        
        # seek the end of file
        f.seek( 0, 2 )
        
        # read the last line and get the length of a line
        while f.read(1) != b'\n':
            
            f.seek( -2, 1 )
            
        lastline = f.readline().decode()
        
        len_line = len( lastline )
        
        lastline = lastline.strip().split()
        
        # if the last row of data is after backup step ?
        # if so, read a previous line; if not, break
        while int(lastline[0]) > backupstep:
            
            f.seek( -2*len_line, 1 )
        
            lastline = f.readline().decode().strip().split() 
        
        # get the backup time for trim ib_force file
        backuptime = float( lastline[1] )
        
        f.truncate()

# --- delete the additional snapshots directories
    
    # get all subfolders in the 'snapshots' directory
    subfolders = [x[0] for x in os.walk( os.getcwd() )]
    subfolders = subfolders[1:]
    
    # select the subfolders after backup step
    subfolders.sort( reverse=True )
    
    for i, subfolder in enumerate( subfolders ):
        
        if int(subfolder[-8:]) <= backupstep:
            break
        
        pass
    
    subfolders = subfolders[0:i]
    
    # if number of subfolders to be deleted > 0, execute cleaning.
    # else raise error.
    
    n_subfolder = len( subfolders )
    
    if n_subfolder == 0:
        
        raise ValueError("There is no snapshots should be deleted!")

    elif n_subfolder > 0:
        
        # delete files first
        for subfolder in subfolders:
            
            os.chdir( subfolder )
            
            for file in os.listdir( subfolder ):
                
                os.remove( file )
        
        # delete subfolder directories        
        for subfolder in subfolders:
            
            os.rmdir( subfolder )

        print( "%d snapshots have been cleared.\n"%n_subfolder )



# ----------------------------------------------------------------------
# >>> Cleaning Probes                                             ( 2 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/12/07  - created
#
# Desc
#
# - walk through every probe file and trim data after a certain step
#
# ----------------------------------------------------------------------

if trim_probe:

    os.chdir( folder )

    os.chdir( 'probes' )

    print( os.getcwd() )

    probefiles = os.listdir( os.getcwd() )

    probefiles.sort()

    # operate each probefiles

    for probefile in probefiles: 
        
        with open( probefile, 'r+b') as f:
            
            # seek the end of file
            f.seek( 0, 2 )
            
            # read the last line and get the length of a line 
            while f.read(1) != b'\n':
                
                f.seek( -2, 1 )
                
            lastline = f.readline().decode()
            
            len_line = len( lastline )

            lastline = lastline.strip().split()
            
            # if the last row of data is after backup step ?
            # if so, read a previous line; if not, break
            # lastline[0] is timestep
            
            while int(lastline[0]) > backupstep:
                
                f.seek( -2*len_line, 1 )
                
                lastline = f.readline().decode().strip().split()

            # truncate the file after backup step
            
            f.truncate()
    
    print( 'probe files have been cleaned.\n' )


# ----------------------------------------------------------------------
# >>> Clean IB Force File                                         ( 3 )
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2022/12/07  - created
#
# Desc
# 
# - trim ib force data after BACKUP TIME
# - backup time is obtained from 'trim_snapshots' part
#
# ----------------------------------------------------------------------

if trim_ibforce:
    
    os.chdir( folder )
    
    print( os.getcwd() )
    
    if 'backuptime' not in globals():
        raise ValueError("snapshot_index.dat file should be included in order to clean ib force files.")
    
    # get all ib force files, only retain filenames instead of paths.
    ibfiles = get_filelist( folder+'/forces', 'force_0' )
    
    for file in ibfiles:
        
        with open( file, 'r+b') as f:
            
            # seek the end of file
            f.seek( 0, 2 )
            
            # read the last line and get the length of a line
            while f.read(1) != b'\n':
                
                f.seek( -2, 1 )
                
            lastline = f.readline().decode()
            
            len_line = len( lastline )
            
            lastline = lastline.strip().split()
            
            # if the last row of data is after backup step ?
            # if so, read a previous line; if not, break
            # lastline[0] is physical time 
            
            while float(lastline[0]) > backuptime:
                
                f.seek( -2*len_line, 1 )
                
                lastline = f.readline().decode().strip().split()

            # truncate the file after backup step
            
            f.truncate()   
            
    print( 'ib force files have been cleaned.' )
