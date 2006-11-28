#!/usr/bin/env python

import sys
import os
import subprocess
import mpi
import Framework.Image
import RO.DS9


# Setup for MPI use
cmd=""
for arg in sys.argv:
    cmd=cmd+arg+" "
# Acquire MPI environmental details
myid_numprocs = mpi.mpi_start(len(sys.argv),cmd)
myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
print "\n-------Start mpiHalfDist_WcsPerMEF on "+str(myid)+"...."

# Acquire test image from local directory or system
try:
    #os.putenv("LSSTProto","/dev/shm")
    #testData = '/dev/shm/tests/WCS/data/642538p.fits'
    testData = os.environ['LSSTProto'] + '/tests/WCS/data/642538p.fits'
except:
    testData = "./642538p.fits"
print "MEF: testData ",testData

# Acquire mosaic configuration file from local directory or system
# 
#    N O T E: both the mosaic and ccd conf files need to be tuned to the
#	      detector specifics
# 
#	      The mosaic conf file defines the CCD conf file to be used.
try:
    mosaicConfFile = os.environ['LSSTProto'] + '/etc/CFHT12K_Mosaic.conf'
    ccdConfFile = os.environ['LSSTProto'] + '/etc/CFHT12K_CCD.conf'
    #mosaicConfFile = '/dev/shm/etc/CFHT12K_Mosaic.conf'
    #ccdConfFile = '/dev/shm/etc/CFHT12K_CCD.conf'
except:
    mosaicConfFile = "./CFHT12K_Mosaic.conf"
    ccdConfFile = "./CFHT12K_CCD.conf"
print "MEF: mosaicConfFile ",mosaicConfFile," ccdConfFile: ",ccdConfFile

# Acquire the Mosaic Image (and relevant HDU info)
immef=Framework.Image.MosaicImage(testData,\
        mosaicConfFile,ccdConfFile)

# determine how many CCD to process
nCCD = immef.NumCCDS()

rootDir = os.getcwd()
print "MEF: Root dir: ",rootDir

# Synchronize before parallel step
mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

# Calculate the WCS for all HDU's at same time
cmd = ("/usr/bin/time mpirun -ssi rpi lamd -D -O -sigs -np "+str(nCCD)+" "+rootDir+"/mpiHalfDist_WcsPerHdu.py "+rootDir)
print "MPI call: "+cmd
res = os.system(cmd)

if not (res == 0) :
    print "Error: Problem building WCS for all "+str(nCCD)+" HDUs in MEF\n"
else:
    print "OK: WCS done for all  "+str(nCCD)+" HDUs in MEF\n"

# Synchronize after step
mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

#-----------------------------------------------------------------------------
# TBD: possibly ..... I don't know what's wanted once the WCS are calculated...
# ?Combine individual FITS files (loaded into subdir) into update MEF?
#
#  ? maybe remove all the transient directories used during calculation?
#-----------------------------------------------------------------------------
     
mpi.mpi_finalize()
print "End mpiHalfDist_WcsPerMEF"

