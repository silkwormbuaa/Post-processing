'''
# ------------------------------------------------------------------------
# Inspect output files from INCA
# ------------------------------------------------------------------------

Description

   This script retrieves the minimum and maximum execution times for all
   tasks recorded in the inca_*.out files within the inca_out directory.

Content

    1. main                     - main routine
    2. assess_performance       - process performance data

History

   Jul 2021  LL  created
   Jan 2024  LL  improved and cleaned

Info

   Refer questions/comments to L. Laguarda at L.LaguardaSanchez@tudelft.nl

'''

# Packages
import struct
import numpy                   as np
import sys
import os

# ------------------------------------------------------------------------
# MAIN                                                               ( 1 )
#
# Author
#  L. Laguarda
#
# History
#  Jun 2021  LL  created
#
# Info
#
# ------------------------------------------------------------------------
def main():

   # Target directory
   out_dir = os.getcwd()

   # Path of target files
   # --------------------

   # List subdirectories
   subdirs = [ out_dir ]
   subdirs.extend( [ f.path for f in os.scandir( out_dir ) if f.is_dir() ] )

   # List inca_*.out files and store their path
   files = []
   for s, subdir in enumerate( subdirs ):

      # List with all files
      files.extend( [ f.path for f in os.scandir( subdir ) if ( f.is_file() and f.name[:5] == 'inca_' ) ] )

   # Inpsect output files
   # --------------------
   content = []
   for i, file in enumerate( files ):
      with open( file ) as f:
         content.extend( [ line.split( '|' ) for line in ( f.readlines() ) if len( line ) > 100 and line[28] == '|' and line[1:5] != 'Name' ] )

   # Reorder list alphabetically
   content.sort()

   # Process data
   assess_performance( content )

# ------------------------------------------------------------------------
# PROCESS PERFORMANCE DATA                                           ( 2 )
#
# Author
#  L. Laguarda
#
# History
#  Jun 2021  LL  created
#
# Info
#
# ------------------------------------------------------------------------
def assess_performance( data ):

   # Retrieve data of interest

   count   = -1
   results = []

   for i, entry in enumerate( data ):

      task = entry[0].strip( ' ' )
      mint = float( entry[2] )
      maxt = float( entry[4] )

      if i == 0 or task != results[count][0]:

         results.append( [ task, mint, maxt, maxt/mint ] ); count += 1

      else:

         results[count][1] = min( results[count][1], mint )
         results[count][2] = max( results[count][2], maxt )
         results[count][3] = results[count][2]/results[count][1]

   # Reorder list alphabetically
   results.sort( key = lambda x: x[3], reverse = True )

   # Output results
   print( '\n\033[92m Performance assessment of target inca_*.out files \033[0m\n' )
   print( ' Task name'.ljust(25) + 't_min [s]'.rjust(15) + 't_max [s]'.rjust(15) + 't_max/t_min'.rjust(15) )
   print( '----------------------------------------------------------------------' )

   for i, result in enumerate( results ):

      # Build output
      string  = ' ' + result[0].ljust(24)
      string += '%15.4E' %( result[1] )
      string += '%15.4E' %( result[2] )
      string += '%15.2f' %( result[3] )
      print( string )
   print( '' )

# ------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------
if __name__ == '__main__': sys.exit( main() )