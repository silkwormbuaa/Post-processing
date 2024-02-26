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
   Feb 2024  WW  assess the average performance of the processors

Info

   Refer questions/comments to L. Laguarda at L.LaguardaSanchez@tudelft.nl

'''

# Packages
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

   count   = -1    # counter for the number of tasks
   results = []

   for i in range( len(data) ):
      
      entry = data[i]
      task  = entry[0].strip( ' ' )
          
      if i == 0 or task != current_task: # data are sorted alphabetically
         
         current_task = task
         freq  = int( entry[1] )
         tmins = [float( entry[2] )]
         tavgs = [float( entry[3] )]
         tmaxs = [float( entry[4] )]
         
      else:

         tmins.append( float( entry[2] ) )
         tavgs.append( float( entry[3] ) )
         tmaxs.append( float( entry[4] ) )
         
      if i+1 > len(data)-1 or data[i+1][0].strip( ' ' ) != current_task:

         tminavg = np.mean(tmins)
         tminr0  = min(tmins) / tminavg
         tminr1  = max(tmins) / tminavg
         
         tmaxavg = np.mean(tmaxs)
         tmaxr0  = min(tmaxs) / tmaxavg
         tmaxr1  = max(tmaxs) / tmaxavg
         
         tavgavg = np.mean(tavgs)
         tavgr0  = min(tavgs) / tavgavg
         tavgr1  = max(tavgs) / tavgavg
         
         results.append( [ task, freq, 
                           tminavg, tminr0, tminr1,
                           tmaxavg, tmaxr0, tmaxr1,
                           tavgavg, tavgr0, tavgr1] )
         
         count += 1


   # Output results
   print( '\n\033[92m Performance assessment of target inca_*.out files \033[0m\n' )
   print( ' Task_name'.ljust(25) + 'freq'.rjust(10), end='' )
   print( 'avg(t_min)'.rjust(15) + 'range_min'.rjust(15) + 'range_max'.rjust(15), end='' )
   print( 'avg(t_max)'.rjust(15) + 'range_min'.rjust(15) + 'range_max'.rjust(15), end='' )
   print( 'avg(t_avg)'.rjust(15) + 'range_min'.rjust(15) + 'range_max'.rjust(15) )
   print( '---------------------------------------------------------------------------------------------------------------------' )

   for i, result in enumerate( results ):

      # Build output
      string  = ' ' + result[0].ljust(24)
      string += '%10d'   %( result[1] )
      for j in range(2,11):
         string += '%15.7f' %( result[j] )
      print( string )
   print( '' )

# ------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------
if __name__ == '__main__': sys.exit( main() )