#!/usr/bin/env python

# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import sys
import lsst.utils
import lsst.afw.display as afwDisplay
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.utils.logging as logUtils
import lsst.pex.config as pexConfig

display = True

verbosity = 4
logUtils.trace_set_at("lsst.ip.diffim", verbosity)

defDataDir = lsst.utils.getPackageDir('afwdata')
imageProcDir = lsst.utils.getPackageDir('ip_diffim')

if len(sys.argv) == 1:
    defTemplatePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2_tmpl.fits")
    defSciencePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2.fits")
    templateMaskedImage = afwImage.MaskedImageF(defTemplatePath)
    scienceMaskedImage = afwImage.MaskedImageF(defSciencePath)
    bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(512, 512))
    templateMaskedImage = afwImage.MaskedImageF(templateMaskedImage, bbox, origin=afwImage.LOCAL)
    scienceMaskedImage = afwImage.MaskedImageF(scienceMaskedImage, bbox, origin=afwImage.LOCAL)

elif len(sys.argv) == 3:
    defTemplatePath = sys.argv[1]
    defSciencePath = sys.argv[2]
    templateMaskedImage = afwImage.MaskedImageF(defTemplatePath)
    scienceMaskedImage = afwImage.MaskedImageF(defSciencePath)
else:
    sys.exit(1)


configAL = ipDiffim.ImagePsfMatchTask.ConfigClass()
configAL.kernel.name = "AL"
subconfigAL = configAL.kernel.active

configDF = ipDiffim.ImagePsfMatchTask.ConfigClass()
configDF.kernel.name = "DF"
subconfigDF = configDF.kernel.active

configDFr = ipDiffim.ImagePsfMatchTask.ConfigClass()
configDFr.kernel.name = "DF"
subconfigDFr = configDFr.kernel.active

subconfigDF.useRegularization = False
subconfigDFr.useRegularization = True

for param in [["spatialKernelOrder", 0],
              ["spatialBgOrder", 0],
              ["usePcaForSpatialKernel", True],
              ["fitForBackground", True]]:
    exec("subconfigAL.%s = %s" % (param[0], param[1]))
    exec("subconfigDF.%s = %s" % (param[0], param[1]))
    exec("subconfigDFr.%s = %s" % (param[0], param[1]))

detConfig = subconfigAL.detectionConfig
kcDetect = ipDiffim.KernelCandidateDetectionF(pexConfig.makePropertySet(detConfig))
kcDetect.apply(templateMaskedImage, scienceMaskedImage)
footprints = kcDetect.getFootprints()

# delta function
psfmatch1 = ipDiffim.ImagePsfMatchTask(config=configDF)
results1 = psfmatch1.subtractMaskedImages(templateMaskedImage, scienceMaskedImage,
                                          candidateList=footprints)
diffim1 = results1.subtractedMaskedImage
spatialKernel1 = results1.psfMatchingKernel
spatialBg1 = results1.backgroundModel
kernelCellSet1 = results1.kernelCellSet

# alard lupton
psfmatch2 = ipDiffim.ImagePsfMatchTask(config=configAL)
results2 = psfmatch2.subtractMaskedImages(templateMaskedImage, scienceMaskedImage,
                                          candidateList=footprints)
diffim2 = results2.subtractedMaskedImage
spatialKernel2 = results2.psfMatchingKernel
spatialBg2 = results2.backgroundModel
kernelCellSet2 = results2.kernelCellSet

# regularized delta function
psfmatch3 = ipDiffim.ImagePsfMatchTask(config=configDFr)
results3 = psfmatch3.subtractMaskedImages(templateMaskedImage, scienceMaskedImage,
                                          candidateList=footprints)
diffim3 = results3.subtractedMaskedImage
spatialKernel3 = results3.psfMatchingKernel
spatialBg3 = results3.backgroundModel
kernelCellSet3 = results3.kernelCellSet


basisList1 = spatialKernel1.getKernelList()
basisList2 = spatialKernel2.getKernelList()
basisList3 = spatialKernel3.getKernelList()

frame = 1
for idx in range(min(5, len(basisList1))):
    kernel = basisList1[idx]
    im = afwImage.ImageD(spatialKernel1.getDimensions())
    ksum = kernel.computeImage(im, False)
    afwDisplay.Display(frame=frame).mtv(im, title="delta function")
    frame += 1

for idx in range(min(5, len(basisList2))):
    kernel = basisList2[idx]
    im = afwImage.ImageD(spatialKernel2.getDimensions())
    ksum = kernel.computeImage(im, False)
    afwDisplay.Display(frame=frame).mtv(im, title="alard lupton")
    frame += 1


for idx in range(min(5, len(basisList3))):
    kernel = basisList3[idx]
    im = afwImage.ImageD(spatialKernel3.getDimensions())
    ksum = kernel.computeImage(im, False)
    afwDisplay.Display(frame=frame).mtv(im, title="regularized delta function")
    frame += 1
