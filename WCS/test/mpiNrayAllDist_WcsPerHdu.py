#!/usr/bin/env python

import sys
import os
import subprocess
import mpi
import Framework.Image
import Framework.Collection
import Pipeline.WCS
import RO.DS9
from numarray import *
import numarray

# Setup for MPI use
cmd=''
for arg in sys.argv:
    cmd=cmd+arg+' '
# Acquire MPI environmental details
myid_numprocs = mpi.mpi_start(len(sys.argv),cmd)
myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
print 'Start mpiWcsPerHdu for ccd:'+str(myid)+'....'

# Acquire test image from local directory or system
try:
    testData = os.environ['LSSTProto'] + '/SampleData/data/642538p.fits'
except:
    testData = './642538p.fits'

# Acquire mosaic configuration file from local directory or system
# 
#    N O T E: both the mosaic and ccd conf files need to be tuned to the
#	      detector specifics
# 
#	      The mosaic conf file defines the CCD conf file to be used.
try:
    mosaicConfFile = os.environ['LSST_POLICY_DIR'] + '/CFHT12K_Mosaic.conf'
    ccdConfFile = os.environ['LSST_POLICY_DIR'] + '/CFHT12K_CCD.conf'
except:
    mosaicConfFile = './CFHT12K_Mosaic.conf'
    ccdConfFile = './CFHT12K_CCD.conf'

# Acquire the Mosaic Image (and relevant HDU info)
immef=Framework.Image.MosaicImage(testData,\
        mosaicConfFile,ccdConfFile)

# Select CCD to process
try:
    ccd=immef.GetCCDImage(int(myid))
except:
    print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,' Cause1: ',sys.exc_value
    mpi.mpi_finalize()
    sys.exit(1)


# Distribute SCAT data to all processors instead of calculating for each 
#
# At the moment, all ccd's have the SkyRegion of the full mosaic
# so just get the star catalog once

# Extract the source catalog stars and then broadcast to all other processors
if myid == 0:  
    # Extract fiducial stars from ascii download of USNO-B named "txtcat"
    #                           and located in $WCS_CAT/txtcat
    try:
        scat=Framework.Collection.StarCollection( \
                                   ccd.GetSkyRegion(),'txtcat', nStarMax=5000)
    except:
        print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,' Cause2: ',sys.exc_value
        mpi.mpi_finalize()
    	sys.exit(1)

    try:
        scat.SortByMag(truncLen=300)
    except:
        print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,' Cause3: ',sys.exc_value
        mpi.mpi_finalize()
    	sys.exit(1)

    # broadcast star count 
    aCount = array(([0]),'i')
    aCount[0] = scat.nStars
    nCount = mpi.mpi_bcast(aCount,1,mpi.MPI_INT,0,mpi.MPI_COMM_WORLD)

    # broadcast scat star quantities 
    id = mpi.mpi_bcast(scat.idArray,scat.nStars,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    ra = mpi.mpi_bcast(scat.raArray,scat.nStars,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    dec = mpi.mpi_bcast(scat.decArray,scat.nStars,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    raPM = mpi.mpi_bcast(scat.raPMArray,scat.nStars,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    decPM = mpi.mpi_bcast(scat.decPMArray,scat.nStars,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    mag1 = mpi.mpi_bcast(scat.mag1Array,scat.nStars,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    mag2 = mpi.mpi_bcast(scat.mag2Array,scat.nStars,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    flux = mpi.mpi_bcast(scat.fluxArray,scat.nStars,mpi.MPI_INT,0,mpi.MPI_COMM_WORLD)
    
else :
    # Receive broadcast of item count
    aCount = array(([0]),'i')
    nCount = mpi.mpi_bcast(aCount,1,mpi.MPI_INT,0,mpi.MPI_COMM_WORLD)
    ncount = nCount[0]

    # Receive broadcast of items
    idArray = array(([0.0]),'f')
    raArray = array(([0.0]),'f')
    decArray = array(([0.0]),'f')
    raPMArray = array(([0.0]),'f')
    decPMArray = array(([0.0]),'f')
    mag1Array = array(([0.0]),'f')
    mag2Array = array(([0.0]),'f')
    fluxArray = array(([0]),'i')

    id  = mpi.mpi_bcast(idArray,ncount,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    ra  = mpi.mpi_bcast(raArray,ncount,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    dec  = mpi.mpi_bcast(decArray,ncount,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    raPM  = mpi.mpi_bcast(raPMArray,ncount,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    decPM  = mpi.mpi_bcast(decPMArray,ncount,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    mag1  = mpi.mpi_bcast(mag1Array,ncount,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    mag2  = mpi.mpi_bcast(mag2Array,ncount,mpi.MPI_DOUBLE,0,mpi.MPI_COMM_WORLD)
    flux  = mpi.mpi_bcast(fluxArray,ncount,mpi.MPI_INT,0,mpi.MPI_COMM_WORLD)
    #print "Id ra dec pra pdec mag1 mag2 flux:"
    #for i in range(ncount):
    #   print myid,i,id[i],ra[i],dec[i],raPM[i],decPM[i],mag1[i],mag2[i],flux[i]

    # Build  'scat' container
    scat = Framework.Collection.StarCollection(ccd.GetSkyRegion(),'handbuilt',ncount,0,id,ra,dec,raPM,decPM,mag1,mag2,flux)
    #for i in range(ncount):
    #   print myid,i,scat.idArray[i],scat.raArray[i],scat.decArray[i],scat.raPMArray[i],scat.decPMArray[i],scat.mag1Array[i],scat.mag2Array[i],scat.fluxArray[i]
     

# position into pre-built HDU-mapped subdirectories
os.chdir(sys.argv[1])
os.chdir('ccd_'+str(myid))

# Select CCD to process
try:
    ccd=immef.GetCCDImage(myid)
except:
    print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,' Cause4: ',sys.exc_value
    mpi.mpi_finalize()
    sys.exit(1)

# Extract sources from that CCD
try:
    sc=ccd.ExtractSources()
except:
    print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,' Cause5: ',sys.exc_value
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
        wcs=Pipeline.WCS.WCS(match)
    except:
        print 'mpiWcsPerHdu: ccd_',myid,': ',sys.exc_type,' Cause6: ',sys.exc_value
        print 'mpiWcsPerHdu: ccd_',myid,': Failure: no WCS constructed'
        mpi.mpi_finalize()
    	sys.exit(1)
except:
    print 'mpiWcsPerHdu: ccd_',myid,' ',sys.exc_type,' Cause7: ',sys.exc_value,''
    print 'mpiWcsPerHdu: ccd_',myid,': Failure: no stars matched, no WCS constructed'
    mpi.mpi_finalize()
    sys.exit(1)
                                                                                
#-----------------------------------------------------------------------------
print 'End mpiWcsPerHdu  for ccd_'+str(myid)
mpi.mpi_finalize()


