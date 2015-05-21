#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2014 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import sys
import numpy as np

import lsst.utils
import lsst.daf.base               as dafBase
import lsst.afw.table              as afwTable
import lsst.afw.image              as afwImage
import lsst.afw.display.ds9        as ds9
import lsst.meas.algorithms        as measAlg
from lsst.meas.algorithms.detection import SourceDetectionTask
from lsst.ip.diffim import DipoleMeasurementTask, DipoleAnalysis

def loadData(imFile=None):
    """Prepare the data we need to run the example"""

    if imFile is None:
        # Load sample input from disk
        afwdataDir = lsst.utils.getPackageDir('afwdata')
        imFile = os.path.join(afwdataDir, "CFHT", "D4", "cal-53535-i-797722_small_1.fits")
    else:
        if not os.path.isfile(imFile):
            print >> sys.stderr, "Input file %s does not exist" % (imFile)
            sys.exit(1)

    exposure = afwImage.ExposureF(imFile)
    psf = measAlg.SingleGaussianPsf(21, 21, 2)
    exposure.setPsf(psf)

    im = exposure.getMaskedImage().getImage()
    im -= np.median(im.getArray())

    # Create the dipole
    offset = 3
    tmpim = im.getArray()[:-offset,:-offset]
    im.getArray()[offset:,offset:] -= tmpim

    return exposure

def run(args):
    exposure = loadData(args.image)
    if args.debug:
        ds9.mtv(exposure, frame=1)

    schema = afwTable.SourceTable.makeMinimalSchema()
        
    # Create the detection task
    config = SourceDetectionTask.ConfigClass()
    config.thresholdPolarity = "both"
    config.background.isNanSafe = True
    config.thresholdValue = 3
    detectionTask = SourceDetectionTask(config=config, schema=schema)
    
    # And the measurement Task
    config = DipoleMeasurementTask.ConfigClass()
    algMetadata = dafBase.PropertyList()
    schema.addField(DipoleMeasurementTask._ClassificationFlag, "F", "probability of being a dipole")
    measureTask = DipoleMeasurementTask(schema, algMetadata, config=config)

    # Create the output table
    tab = afwTable.SourceTable.make(schema)
    
    # Process the data
    results = detectionTask.run(tab, exposure)

    # Merge the positve and negative sources
    fpSet = results.fpSets.positive
    growFootprint = 2
    fpSet.merge(results.fpSets.negative, growFootprint, growFootprint, False)
    diaSources = afwTable.SourceCatalog(tab)
    fpSet.makeSources(diaSources)

    print "Merged %s Sources into %d diaSources (from %d +ve, %d -ve)" % (len(results.sources), 
        len(diaSources), results.fpSets.numPos, results.fpSets.numNeg)

    measureTask.measure(exposure, diaSources)

    # Display dipoles if debug enabled
    if args.debug:
        dpa = DipoleAnalysis()
        dpa.displayDipoles(exposure, diaSources)
 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Demonstrate the use of SourceDetectionTask and DipoleMeasurement}Task")

    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument("--image", "-i", help="User defined image", default=None)

    args = parser.parse_args()

    if args.debug:
        try:
            import debug
            debug.lsstDebug.frame = 2
        except ImportError as e:
            print >> sys.stderr, e

    run(args)
