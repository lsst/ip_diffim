#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.

from numpy import (random as np_random,
                   array as np_array)

## LSST imports:
import lsst.afw.image as afw_image
import lsst.utils.tests as lsst_tests
from lsst.afw.table import (SourceTable, SourceCatalog)
from lsst.ip.diffim import dipoleFitTask as dft
from lsst.meas.base import SingleFrameMeasurementConfig
from lsst.meas.algorithms import (SourceDetectionConfig, SourceDetectionTask)

posImage = afw_image.ExposureF('calexp-0289820_11.fits')
negImage = afw_image.ExposureF('matchexp-11.fits')
diffim = afw_image.ExposureF('diffexp-11.fits')

def runDetection(exposure):
    # Customize the detection task a bit (optional)
    detectConfig = SourceDetectionConfig() 
    detectConfig.returnOriginalFootprints = False # should be the default 
    detectConfig.thresholdValue = 5.
    detectConfig.reEstimateBackground = True
    detectConfig.thresholdType = "pixel_stdev"

    # Create the detection task. We pass the schema so the task can declare a few flag fields
    schema = SourceTable.makeMinimalSchema()
    detectTask = SourceDetectionTask(config=detectConfig, schema=schema)

    table = SourceTable.make(schema)
    detectResult = detectTask.run(table, exposure)
    posCatalog = detectResult.sources
    return posCatalog

posCatalog = runDetection(posImage)

negImage.setPsf(posImage.getPsf())
negCatalog = runDetection(negImage)

import imp, os
print os.getenv('IP_DIFFIM_DIR')+'/tests/testDipoleFitter.py'
dtUtils = imp.load_source('dtUtils', os.getenv('IP_DIFFIM_DIR')+'/tests/testDipoleFitter.py')
from lsst.ip.diffim import dipoleFitTask as dft

catalog = dtUtils.DipoleTestUtils.detectDipoleSources(
    diffim, posImage, posCatalog, negImage, negCatalog)

for i,s in enumerate(catalog):
    fp = s.getFootprint()
    if (len(fp.getPeaks()) <= 1): continue
    print i, fp.getBBox(), fp.getNpix(), len(fp.getPeaks())
    for pk in fp.getPeaks():
        print '   FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

    result = dft.DipoleFitAlgorithm.fitDipole_new(
        diffim, s, posImage, negImage, rel_weight=0.5, separateNegParams=False,
        verbose=False, display=False)


s = catalog[115]
fp = s.getFootprint()
print fp.getBBox(), fp.getNpix()
for pk in fp.getPeaks():
    print '   FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

#dft.DipoleUtils.displayCutouts(catalog[0], diffim, posImage, negImage)
dft.DipoleFitAlgorithm.fitDipole_new(
    diffim, s, posImage, negImage, rel_weight=1., separateNegParams=False,
    verbose=False, display=True)

## ===============================

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

bbox = afwGeom.Box2I(afwGeom.Point2I(230, 1080), afwGeom.Extent2I(100, 100))
subim = afwImage.ImageF(diffim.getMaskedImage().getImage(), bbox, afwImage.PARENT)
