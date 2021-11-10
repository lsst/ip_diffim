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
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.log.utils as logUtils

afwDisplay.setDefaultMaskTransparency(75)

subBackground = True
display = True
fwhm = 6.8
warp = True

verbosity = 4
logUtils.traceSetAt("lsst.ip.diffim", verbosity)

defDataDir = lsst.utils.getPackageDir('afwdata')
imageProcDir = lsst.utils.getPackageDir('ip_diffim')

if len(sys.argv) == 1:
    defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                  "v26-e0-c011-a10.sci.fits")
    defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                   "v5-e0-c011-a10.sci.fits")
elif len(sys.argv) == 3:
    defSciencePath = sys.argv[1]
    defTemplatePath = sys.argv[2]
else:
    sys.exit(1)

defOutputPath = "diffImage"
templateExposure = afwImage.ExposureF(defTemplatePath)
scienceExposure = afwImage.ExposureF(defSciencePath)

config = ipDiffim.ImagePsfMatchTask.ConfigClass()
config.kernel.name = "AL"
subconfig = config.kernel.active


if warp:
    warper = afwMath.Warper.fromConfig(subconfig.warpingConfig)
    templateExposure = warper.warpExposure(scienceExposure.getWcs(), templateExposure,
                                           destBBox=scienceExposure.getBBox())
if subBackground:
    # Do in AFW
    diffimTools.backgroundSubtract(subconfig.afwBackgroundConfig,
                                   [templateExposure.getMaskedImage(),
                                    scienceExposure.getMaskedImage()])
    subconfig.fitForBackground = False
else:
    # Do in IP_DIFFIM
    subconfig.fitForBackground = True


frame = 0
afwDisplay.Display(frame=frame).mtv(templateExposure, title="Template")
frame += 1
afwDisplay.Display(frame=frame).mtv(scienceExposure, title="Sci Im")

psfmatch = ipDiffim.ImagePsfMatchTask(config=config)
try:
    candidateList = psfmatch.makeCandidateList(templateExposure, scienceExposure, subconfig.kernelSize)
    results = psfmatch.subtractMaskedImages(templateExposure.getMaskedImage(),
                                            scienceExposure.getMaskedImage(), candidateList)
    diffim = results.subtractedImage
    spatialKernel = results.psfMatchingKernel
    spatialBg = results.backgroundModel
    kernelCellSet = results.kernelCellSet
except Exception:
    print('FAIL')
    sys.exit(1)

if False:
    print(spatialKernel.getSpatialFunctionList()[0].toString())
    print(spatialKernel.getKernelParameters())
    print(spatialKernel.getSpatialParameters())
    import pdb
    pdb.set_trace()

# Lets see what we got
if display:
    mos = afwDisplay.utils.Mosaic()

    # Inputs
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False):  # False = include bad candidates
            rchi2 = cand.getChi2()

            # No kernels made
            if cand.getStatus() == afwMath.SpatialCellCandidate.UNKNOWN:
                continue

            try:
                im = cand.getKernelImage(ipDiffim.KernelCandidateF.RECENT)
            except Exception:
                continue

            if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                statStr = "Good"
            elif cand.getStatus() == afwMath.SpatialCellCandidate.BAD:
                statStr = "Bad"
            mos.append(im, "#%d: %.1f (%s)" % (cand.getId(), rchi2, statStr))

    mosaic = mos.makeMosaic()
    frame += 1
    afwDisplay.Display(frame=frame).mtv(mosaic, title="Raw Kernels")
    mos.drawLabels(frame=frame)

    # KernelCandidates
    frame += 1
    diffimTools.displayCandidateResults(kernelCellSet, frame)

    # Bases
    mos.reset()
    basisList = spatialKernel.getKernelList()
    for idx in range(len(basisList)):
        kernel = basisList[idx]
        im = afwImage.ImageD(spatialKernel.getDimensions())
        ksum = kernel.computeImage(im, False)
        mos.append(im, "K%d" % (idx))
    mosaic = mos.makeMosaic()
    frame += 1
    afwDisplay.Display(frame=frame).mtv(mosaic, title="Kernel Bases")
    mos.drawLabels(frame=frame)

    # Spatial model
    mos.reset()
    width = templateExposure.getWidth()
    height = templateExposure.getHeight()
    stamps = []
    stampInfo = []
    for x in (0, width//2, width):
        for y in (0, height//2, height):
            im = afwImage.ImageD(spatialKernel.getDimensions())
            ksum = spatialKernel.computeImage(im,
                                              False,
                                              afwImage.indexToPosition(x),
                                              afwImage.indexToPosition(y))
            mos.append(im, "x=%d y=%d kSum=%.2f" % (x, y, ksum))

    mosaic = mos.makeMosaic()
    frame += 1
    afwDisplay.Display(frame=frame).mtv(mosaic, title="Spatial Kernels")
    mos.drawLabels(frame=frame)

    # Background
    backgroundIm = afwImage.ImageF(afwGeom.Extent2I(templateExposure.getWidth(),
                                   templateExposure.getHeight()), 0)
    backgroundIm += spatialBg
    frame += 1
    afwDisplay.Display(frame=frame).mtv(backgroundIm, title="Background model")

    # Diffim!
    diffIm = ipDiffim.convolveAndSubtract(templateExposure.getMaskedImage(),
                                          scienceExposure.getMaskedImage(),
                                          spatialKernel,
                                          spatialBg)
    frame += 1
    afwDisplay.Display(frame=frame).mtv(diffIm, title="Diffim")


# examples/runSpatialModel.py $AFWDATA_DIR/DC3a-Sim/sci/v5-e0/v5-e0-c011-a00.sci
# ... $AFWDATA_DIR/DC3a-Sim/sci/v26-e0/v26-e0-c011-a00.sci
