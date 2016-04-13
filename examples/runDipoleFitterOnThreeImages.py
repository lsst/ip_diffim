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

import numpy as np
from matplotlib import pyplot as plt

#  LSST imports:
import lsst.afw.image as afw_image
from lsst.afw.table import (SourceTable, SourceCatalog)
from lsst.meas.base import SingleFrameMeasurementConfig
from lsst.meas.algorithms import (SourceDetectionConfig, SourceDetectionTask)

posImage = afw_image.ExposureF('calexp-0289820_11.fits')
negImage = afw_image.ExposureF('matchexp-11.fits')
diffim = afw_image.ExposureF('diffexp-11.fits')

negImage.setPsf(posImage.getPsf())

import imp, os
print os.getenv('IP_DIFFIM_DIR')+'/tests/testDipoleFitter.py'
dtUtils = imp.load_source('dtUtils', os.getenv('IP_DIFFIM_DIR')+'/tests/testDipoleFitter.py')
from lsst.ip.diffim import dipoleFitTask as dft

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
negCatalog = runDetection(negImage)

#### Now run the DipoleFitTask plugin on the image.

from lsst.afw.table import (SourceTable, SourceCatalog)
from lsst.meas.algorithms import (SourceDetectionConfig, SourceDetectionTask)

## Create the various tasks and schema -- avoid code reuse.
img = dtUtils.DipoleTestImage()
img.diffim = diffim
img.posImage = posImage
img.negImage = negImage
detectTask, schema = img.detectDipoleSources(detectSigma=5.5, doMerge=False)

measureConfig = SingleFrameMeasurementConfig()

# Modify the set of active plugins ('.names' behaves like a Python set)
#measureConfig.plugins.names.remove("base_PsfFlux")
#measureConfig.plugins.names.remove("base_GaussianCentroid")
#measureConfig.plugins.names.remove("base_SdssCentroid")
#measureConfig.plugins.names.remove("base_GaussianFlux")
#measureConfig.plugins.names.remove("base_SdssShape")

measureConfig.slots.modelFlux = "ip_diffim_DipoleFit"

measureConfig.plugins.names |= ["base_CircularApertureFlux",
                                "base_PixelFlags",
                                "base_SkyCoord",
                                "base_PsfFlux",
                                "ip_diffim_NaiveDipoleCentroid",
                                "ip_diffim_NaiveDipoleFlux",
                                "ip_diffim_PsfDipoleFlux"]

# Disable aperture correction, which requires having an ApCorrMap attached to
# the Exposure (it'll warn if it's not present and we don't explicitly disable it).
measureConfig.doApplyApCorr = "no"

# Here is where we make the dipole fitting task. It can run the other measurements as well.
measureTask = dft.DipoleFitTask(config=measureConfig, schema=schema)

# Erase existing detection mask planes
mask  = diffim.getMaskedImage().getMask()
mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

table = SourceTable.make(schema)
detectResult = detectTask.run(table, diffim)
catalog = detectResult.sources

fpSet = detectResult.fpSets.positive
grow = 5
fpSet.merge(detectResult.fpSets.negative, grow, grow, False)
sources = SourceCatalog(table)
fpSet.makeSources(sources)

for i,s in enumerate(sources):
    fp = s.getFootprint()
    print i, s.getId(), fp.getBBox(), fp.getNpix(), len(fp.getPeaks())
    for pk in fp.getPeaks():
        print '   FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

measureTask.run(sources, diffim, posImage, negImage)


fitResults = []
for i,s in enumerate(sources):
    fp = s.getFootprint()
    if (len(fp.getPeaks()) <= 1): continue
    print i, s.getId(), fp.getBBox(), fp.getNpix(), len(fp.getPeaks())
    print '   ', s.extract("ip_diffim_DipoleFit_flag_classification")
    for pk in fp.getPeaks():
        print '   FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()
    fitres = img.fitDipoleSource(s, verbose=False, display=False)
    fitResults.append(fitres)


s = dft.DipolePlotUtils.searchCatalog(sources, 617, 209)
print s.extract("ip_diffim_DipoleFit*")

if False:
    img.displayCutouts(s, False)
    fitResult = img.fitDipoleSource(s, verbose=True, display=True)

#%timeit img.fitDipoleSource(s, verbose=False, display=False)
fitResult = img.fitDipoleSource(s, verbose=True, display=False)

if False:
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('multipage_pdf.pdf')
    page = 1
    for i,s in enumerate(sources):
        isDipole = s.extract("ip_diffim_DipoleFit*")['ip_diffim_DipoleFit_flag_classification']
        if isDipole or s.getId() == 331:
            print i, s.getId(), page, isDipole
            page += 1
            fig = dft.DipolePlotUtils.displayCutouts(s, diffim, posImage, negImage,
                                                     asHeavyFootprint=True, title='%d %d'%(i, s.getId()))
            pdf.savefig(fig)
    pdf.close()

if False:
    ##s = dft.DipolePlotUtils.searchCatalog(sources, 1592,900)
    s = dft.DipolePlotUtils.searchCatalog(sources, 1432, 1077)
    result = dft.DipoleFitAlgorithm.fitDipole_new(
        diffim, s, posImage, negImage, rel_weight=0.5, separateNegParams=False,
        verbose=True, display=True)
    print result.psfFitOrientation, \
        np.sqrt((result.psfFitPosCentroidX - result.psfFitNegCentroidX)**2. + \
                (result.psfFitPosCentroidY - result.psfFitNegCentroidY)**2.)
