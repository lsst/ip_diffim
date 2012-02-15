#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
import eups
import lsst.afw.geom as afwGeom
import lsst.afw.image.imageLib as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLogging
import lsst.pex.config as pexConfig
import lsst.afw.display.ds9 as ds9

display = True

verbosity = 5
pexLogging.Trace_setVerbosity("lsst.ip.diffim", verbosity)

defDataDir   = eups.productDir("afwdata") 
imageProcDir = eups.productDir("ip_diffim")

if len(sys.argv) == 1:
    defTemplatePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2_tmpl")
    defSciencePath  = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2")
    templateMaskedImage = afwImage.MaskedImageF(defTemplatePath)
    scienceMaskedImage  = afwImage.MaskedImageF(defSciencePath)
    bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(512, 512))
    templateMaskedImage = afwImage.MaskedImageF(templateMaskedImage, bbox, afwImage.LOCAL)
    scienceMaskedImage  = afwImage.MaskedImageF(scienceMaskedImage, bbox, afwImage.LOCAL)
    
elif len(sys.argv) == 3:
    defTemplatePath = sys.argv[1]
    defSciencePath  = sys.argv[2]
    templateMaskedImage = afwImage.MaskedImageF(defTemplatePath)
    scienceMaskedImage  = afwImage.MaskedImageF(defSciencePath)
else:
    sys.exit(1)
    

configAL  = ipDiffim.PsfMatchConfigAL()
configDF  = ipDiffim.PsfMatchConfigDF()
configDFr = ipDiffim.PsfMatchConfigDF()

configDF.useRegularization  = False
configDFr.useRegularization = True

for param in [["spatialKernelOrder", 0],
              ["spatialBgOrder", 0],
              ["usePcaForSpatialKernel", True],
              ["fitForBackground", True]]:
    exec("configAL.%s = %s" % (param[0], param[1]))
    exec("configDF.%s = %s" % (param[0], param[1]))
    exec("configDFr.%s = %s" % (param[0], param[1]))

detConfig = configAL.detectionConfig
kcDetect = ipDiffim.KernelCandidateDetectionF(pexConfig.makePolicy(detConfig))
kcDetect.apply(templateMaskedImage, scienceMaskedImage)
footprints = kcDetect.getFootprints()

# delta function
psfmatch1 = ipDiffim.ImagePsfMatch(configDF)
results1  = psfmatch1.subtractMaskedImages(templateMaskedImage, scienceMaskedImage, footprints = footprints)
diffim1, spatialKernel1, spatialBg1, kernelCellSet1 = results1

# alard lupton
psfmatch2 = ipDiffim.ImagePsfMatch(configAL)
results2  = psfmatch2.subtractMaskedImages(templateMaskedImage, scienceMaskedImage, footprints = footprints)
diffim2, spatialKernel2, spatialBg2, kernelCellSet2 = results2

# regularized delta function
psfmatch3 = ipDiffim.ImagePsfMatch(configDFr)
results3  = psfmatch3.subtractMaskedImages(templateMaskedImage, scienceMaskedImage, footprints = footprints)
diffim3, spatialKernel3, spatialBg3, kernelCellSet3 = results3


basisList1 = spatialKernel1.getKernelList()
basisList2 = spatialKernel2.getKernelList()
basisList3 = spatialKernel3.getKernelList()

frame = 1
for idx in range(min(5, len(basisList1))):
    kernel = basisList1[idx]
    im     = afwImage.ImageD(spatialKernel1.getDimensions())
    ksum   = kernel.computeImage(im, False)    
    ds9.mtv(im, frame=frame)
    frame += 1

for idx in range(min(5, len(basisList2))):
    kernel = basisList2[idx]
    im     = afwImage.ImageD(spatialKernel2.getDimensions())
    ksum   = kernel.computeImage(im, False)    
    ds9.mtv(im, frame=frame)
    frame += 1


for idx in range(min(5, len(basisList3))):
    kernel = basisList3[idx]
    im     = afwImage.ImageD(spatialKernel3.getDimensions())
    ksum   = kernel.computeImage(im, False)    
    ds9.mtv(im, frame=frame)
    frame += 1

