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

## Here is an example of how to run the algorithm on each source.
## Below we will just use the DipoleFitTask

for i,s in enumerate(catalog):
    fp = s.getFootprint()
    if (len(fp.getPeaks()) <= 1): continue
    print i, fp.getBBox(), fp.getNpix(), len(fp.getPeaks())
    for pk in fp.getPeaks():
        print '   FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

    result = dft.DipoleFitAlgorithm.fitDipole_new(
        diffim, s, posImage, negImage, rel_weight=0.5, separateNegParams=False,
        verbose=False, display=False)

from matplotlib import pyplot as plt

s = catalog[176]
fp = s.getFootprint()
print fp.getBBox(), fp.getNpix()
for pk in fp.getPeaks():
    print '   FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

dft.DipoleFitAlgorithm.fitDipole_new(
    diffim, s, posImage, negImage, rel_weight=0.5, separateNegParams=False,
    verbose=True, display=True)

dft.DipoleFitAlgorithm.fitDipole_new(
    diffim, s, rel_weight=0., separateNegParams=False,
    verbose=True, display=True)

#catalog = dft.DipolePlotUtils.makeHeavyCatalog(catalog, diffim)

dft.DipolePlotUtils.displayCutouts(s, diffim, posImage, negImage)
plt.show()
dft.DipolePlotUtils.displayCutouts(s, diffim, posImage, negImage, asHeavyFootprint=True)
plt.show()

s = dft.DipolePlotUtils.searchCatalog(catalog, 1354, 825)
dft.DipolePlotUtils.displayCutouts(s, diffim, posImage, negImage, asHeavyFootprint=True)
plt.show()




#### Now run the DipoleFitTask plugin on the image.

from lsst.afw.table import (SourceTable, SourceCatalog)
from lsst.meas.algorithms import (SourceDetectionConfig, SourceDetectionTask)

## Create the various tasks and schema -- avoid code reuse.
detectTask, deblendTask, schema = dtUtils.DipoleTestUtils.detectDipoleSources(
    diffim, posImage, posCatalog, negImage, negCatalog, doMerge=False)

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

table = SourceTable.make(schema)
detectResult = detectTask.run(table, diffim)
catalog = detectResult.sources
deblendTask.run(diffim, catalog, psf=diffim.getPsf())

fpSet = detectResult.fpSets.positive
fpSet.merge(detectResult.fpSets.negative, 2, 2, False)
sources = SourceCatalog(table)
fpSet.makeSources(sources)

for i,s in enumerate(sources):
    fp = s.getFootprint()
    if (len(fp.getPeaks()) <= 1): continue
    print i, fp.getBBox(), fp.getNpix(), len(fp.getPeaks())
    for pk in fp.getPeaks():
        print '   FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

measureTask.run(sources, diffim, posImage, negImage)

s = sources[0]
print s.extract("ip_diffim_DipoleFit*")

if False:
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('multipage_pdf.pdf')
    page = 1
    for i,s in enumerate(sources):
        isDipole = s.extract("ip_diffim_DipoleFit*")['ip_diffim_DipoleFit_flag_classification']
        if isDipole:
            print i, page, isDipole
            page += 1
            fig = dft.DipolePlotUtils.displayCutouts(s, diffim, posImage, negImage, asHeavyFootprint=True)
            pdf.savefig(fig)
    pdf.close()
