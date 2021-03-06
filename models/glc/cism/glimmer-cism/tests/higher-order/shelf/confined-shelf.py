#!/usr/bin/env python
# This script runs a "Confined Shelf Experiment".
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Written by Glen Granzow at the University of Montana on April 8, 2010

# Slight alterations by SFP on 2-3-11 for Glimmer-CISM 2.0 relase

import sys, os, glob, shutil, numpy
from netCDF import *
from ConfigParser import ConfigParser

# Parse command-line options
from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-c", "--config", dest="configfile", type='string', default='confined-shelf.config', help="Name of .config file to use for the run", metavar="FILE")
optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel [default: perform a serial run]', metavar="NUMPROCS")
optparser.add_option('-e','--exec',dest='executable',default='./cism_driver',help='Set path to the CISM executable')
for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()


# Create a netCDF file according to the information in the config file.
parser = ConfigParser()
parser.read(options.configfile)
nx = int(parser.get('grid','ewn'))
ny = int(parser.get('grid','nsn'))
nz = int(parser.get('grid','upn'))
dx = float(parser.get('grid','dew'))
dy = float(parser.get('grid','dns'))
filename = parser.get('CF input', 'name')

print 'Writing', filename
try:
  netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
except TypeError:
  netCDFfile = NetCDFFile(filename,'w')

netCDFfile.createDimension('time',1)
netCDFfile.createDimension('x1',nx)
netCDFfile.createDimension('y1',ny)
netCDFfile.createDimension('level',nz)
netCDFfile.createDimension('x0',nx-1) # staggered grid 
netCDFfile.createDimension('y0',ny-1)

x = dx*numpy.arange(nx,dtype='float32')
y = dx*numpy.arange(ny,dtype='float32')

netCDFfile.createVariable('time','f',('time',))[:] = [0]
netCDFfile.createVariable('x1','f',('x1',))[:] = x.tolist()
netCDFfile.createVariable('y1','f',('y1',))[:] = y.tolist()

netCDFfile.createVariable('x0','f',('x0',))[:] = (dx/2 + x[:-1]).tolist()
netCDFfile.createVariable('y0','f',('y0',))[:] = (dy/2 + y[:-1]).tolist()

# *SFP* this has been changed so that the default value for 'flwa' is the same 
# as in the EISMINT-shelf test documentation, tests 3 & 4, found at:
# http://homepages.vub.ac.be/~phuybrec/eismint/iceshelf.html

## Check to make sure that the flow law parameter in the config file is correct.
#default_flwa = float(parser.get('parameters','default_flwa'))
#if default_flwa != 4.6e-18:
#  print 'WARNING: The parameter default_flwa in',options.configfile,'should be 4.6e-18'
#  print '         Currently it is',default_flwa

# *SFP* removed periodic option
# Determine from the config file whether periodic boundary conditions are to be
# imposed in the x direction.  
#periodic_ew = int(parser.get('options','periodic_ew'))

# Calculate values for the required variables.
thk  = numpy.zeros([1,ny,nx],dtype='float32')
topg  = numpy.zeros([1,ny,nx],dtype='float32')
beta = numpy.empty([1,ny-1,nx-1],dtype='float32')
kbc  = numpy.zeros([1,ny-1,nx-1],dtype='int')
acab = numpy.zeros([1,ny,nx],dtype='float32') # *sfp* added acab field for prog. runs 
zero = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

# *SFP* added topg var so that areas of no slip are consistent w/ grounded ice on bedrock
topg[:] = -2000.0

# *SFP* changed to be in line w/ EISMINT-shelf tests 3&4 

# shelf bc applied at bottom (DEFAULT FOR TEST CASE - other options below for testing bcs)
thk[0,4:-2,2:-2] = 500.     
kbc[0,ny-4:,:]  = 1
kbc[0,:,:3] = 1
kbc[0,:,nx-4:] = 1
topg[0,ny-4:,:]  = -440 
topg[0,:,:4] = -440
topg[0,:,nx-4:] = -440

# shelf bc applied at top    
#thk[0,2:-4,2:-2] = 500.     
#kbc[0,:3,:]  = 1
#kbc[0,:,:3] = 1
#kbc[0,:,nx-4:] = 1
#topg[0,:4,:]  = -440
#topg[0,:,:4] = -440
#topg[0,:,nx-4:] = -440

# shelf bc applied at right     ! NOTE that shelf is wider slightly wider in ns than in ew direction  
#thk[0,2:-2,2:-4] = 500.     
#kbc[0,:,:3]  = 1
#kbc[0,:3,:] = 1
#kbc[0,ny-4:,:] = 1
#topg[0,:,:4]  = -440
#topg[0,:4,:] = -440
#topg[0,ny-4:,:] = -440

# shelf bc applied at left     ! NOTE that shelf is wider slightly wider in ns than in ew direction  
#thk[0,2:-2,4:-2] = 500.     
#kbc[0,:,nx-4:]  = 1
#kbc[0,:3,:] = 1
#kbc[0,ny-4:,:] = 1
#topg[0,:,nx-4:]  = -440
#topg[0,:4,:] = -440
#topg[0,ny-4:,:] = -440

#if not periodic_ew:    *SFP* removed periodic option

beta[0,:,:] = 0 

acab[:] = 0.25
acab[0,ny-3:,:]  = 0    # zero out accum at edges to avoid buildup where u=0
acab[0,:,:3] = 0
acab[0,:,nx-3:] = 0

# *SFP* calculate stream profile for upstream end
for i in range(nx-2):
  x = float( i ) / (nx-2) - 0.5
  vvel[0,:,ny-4,i] = -1.5e3 * 1/(2*3.141592654*0.125) * numpy.exp( -x**2 / (2*0.125**2) )

# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk',      'f',('time','y1','x1'))[:] = thk.tolist()
netCDFfile.createVariable('acab',     'f',('time','y1','x1'))[:] = acab.tolist()
netCDFfile.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kbc.tolist()
netCDFfile.createVariable('topg',     'f',('time','y1','x1'))[:] = topg.tolist()
netCDFfile.createVariable('beta',     'f',('time','y0','x0'))[:] = beta.tolist()
netCDFfile.createVariable('uvel',  'f',('time','level','y0','x0'))[:] = zero.tolist()

# *sfp* first option below adds ice stream vel profile for kin bc at upstream end
# *sfp* ... comment out for standard test case
#netCDFfile.createVariable('vvel',  'f',('time','level','y0','x0'))[:] = vvel.tolist()
netCDFfile.createVariable('vvel',  'f',('time','level','y0','x0'))[:] = zero.tolist()

netCDFfile.close()


# =====================================
# Run CISM
print 'Running CISM for the confined-shelf experiment'
print '==============================================\n'
if options.parallel == None:
   # Perform a serial run
   runstring = options.executable + ' ' + options.configfile
   print 'Executing serial run with:  ' + runstring + '\n\n'
   os.system(runstring)
else:
   # Perform a parallel run
   if options.parallel <= 0:
      sys.exit( 'Error: Number of processors specified for parallel run is <=0.' )
   else:
      # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure
      if os.system('which openmpirun > /dev/null') == 0:
         mpiexec = 'openmpirun -np ' + str(options.parallel)
      elif os.system('which mpirun > /dev/null') == 0:
         mpiexec = 'mpirun -np ' + str(options.parallel)
      elif os.system('which aprun > /dev/null') == 0:
         mpiexec = 'aprun -n ' + str(options.parallel)
      elif os.system('which mpirun.lsf > /dev/null') == 0:
         # mpirun.lsf does NOT need the number of processors (options.parallel)
         mpiexec = 'mpirun.lsf'
      else:
         sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./cism_driver confined-shelf.config')
      runstring = mpiexec + ' ' + options.executable + ' ' + options.configfile
      print 'Executing parallel run with:  ' + runstring + '\n\n'
      os.system(runstring)  # Here is where the parallel run is actually executed!


