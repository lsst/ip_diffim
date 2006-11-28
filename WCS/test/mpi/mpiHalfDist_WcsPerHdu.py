#!/usr/bin/env python

import sys
import os
import subprocess
import mpi
import Framework.Image
import Framework.Collection
import Pipeline.WCS
import RO.DS9

# Setup for MPI use
cmd=''
for arg in sys.argv:
    cmd=cmd+arg+' '
# Acquire MPI environmental details
myid_numprocs = mpi.mpi_start(len(sys.argv),cmd)
myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
print '\nStart mpiWcsPerHdu for ccd:'+str(myid)+'....'

# Acquire test image from local directory or system
try:
    #os.putenv("LSSTProto","/dev/shm")
    #testData = '/dev/shm/tests/WCS/data/642538p.fits'
    testData = os.environ['LSSTProto'] + '/tests/WCS/data/642538p.fits'
except:
    testData = './642538p.fits'
print "Hdu: testData ",testData


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
    mosaicConfFile = './CFHT12K_Mosaic.conf'
    ccdConfFile = './CFHT12K_CCD.conf'
print "Hdu: mosaicConfFile ",mosaicConfFile," ccdConfFile: ",ccdConfFile

# Acquire the Mosaic Image (and relevant HDU info)
immef=Framework.Image.MosaicImage(testData,\
        mosaicConfFile,ccdConfFile)

# Select CCD to process
try:
    ccd=immef.GetCCDImage(int(myid))
except:
    print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,'\n Cause1: ',sys.exc_value,'\n'
    mpi.mpi_finalize()
    sys.exit(1)

# Build subdir-per-CCD to ensure no file i/o collisions due to parallel proc
rootDir = sys.argv[1]
os.chdir(sys.argv[1])
print "Hdu: Root dir: ",rootDir
ccdDir = "ccd_"+str(myid)
try:
    os.mkdir(ccdDir)
except OSError:
    pass

# position into pre-built HDU-mapped subdirectories
os.chdir('ccd_'+str(myid))

## ------ TBD: pass the SCAT data around instead of calc for each CCD ------
## At the moment, all ccd's have the SkyRegion of the full mosaic
## so just get the star catalog once
if 0 == 0:  #Change to myid==0 when computing once and sharing results
    # Extract fiducial stars from ascii download of USNO-B named "txtcat"
    #                           and located in $WCS_CAT/txtcat
    try:
        scat=Framework.Collection.StarCollection( \
                                   ccd.GetSkyRegion(),'txtcat', nStarMax=5000)
    except:
        print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,'\nCause2: ',sys.exc_value,'\n'
        mpi.mpi_finalize()
    	sys.exit(1)

    try:
        scat.SortByMag(truncLen=300)
    except:
        print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,'\nCause3: ',sys.exc_value,'\n'
        mpi.mpi_finalize()
    	sys.exit(1)

    #print "scat Id ra dec pra pdec mag1 mag2 flux:"
    for i in range(scat.nStars):
    	print myid,i,scat.idArray[i],scat.raArray[i],scat.decArray[i],scat.raPMArray[i],scat.decPMArray[i],scat.mag1Array[i],scat.mag2Array[i],scat.fluxArray[i]

# TBD
#    # mpi send the scat to all processors
#else:
#    # mpi receive the scat info



# Select CCD to process
try:
    ccd=immef.GetCCDImage(myid)
except:
    print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,'\n Cause4: ',sys.exc_value,'\n'
    mpi.mpi_finalize()
    sys.exit(1)

# Extract sources from that CCD
try:
    sc=ccd.ExtractSources()
except:
    print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,'\n Cause5: ',sys.exc_value,'\n'
    mpi.mpi_finalize()
    sys.exit(1)

                                                                                
# now do WCS processing extracted from mpiWcsPerHdu.py

# Initialize the match class with the source and fiducial star collections
match=Pipeline.WCS.StarSourceMatchCollection(scat, sc, ccd.GetMetaData())
                                                                                
# Perform the star match and acquire the transformation parameters
# match the source and fiducial stars; results remain within match instance
try:
    match.StarMatch()
    # Create the WCS for the CCD
    try:
        wcs=Pipeline.WCS(match)
	#print 'mpiWcsPerHdu: ccd_',myid,': world coords:',wcs.worldCoor,'\n'
    except:
        print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,'\nCause6: ',sys.exc_value,'\n'
        print 'mpiWcsPerHdu: ccd_',myid,': Failure: no WCS constructed\n'
        mpi.mpi_finalize()
    	sys.exit(1)
except:
    print 'mpiWcsPerHdu: ccd_',myid,' ',sys.exc_type,'\nCause7: ',sys.exc_value,'\n'
    print 'mpiWcsPerHdu: ccd_',myid,': Failure: no stars matched, no WCS constructed\n'
    mpi.mpi_finalize()
    sys.exit(1)
                                                                                
#-----------------------------------------------------------------------------
print 'End mpiWcsPerHdu  for ccd_'+str(myid)
mpi.mpi_finalize()


