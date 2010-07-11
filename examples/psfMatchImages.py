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

import sys, os
import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLog

verbosity = 3
pexLog.Trace_setVerbosity('lsst.ip.diffim', verbosity)

imageProcDir  = eups.productDir("ip_diffim")
defPolicyPath = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")

imageToConvolve     = afwImage.MaskedImageF(sys.argv[1])
imageToNotConvolve  = afwImage.MaskedImageF(sys.argv[2])
outputImage         = sys.argv[3]

policy              = ipDiffim.generateDefaultPolicy(defPolicyPath)
policy.set("detThreshold", 5.0)
policy.set("kernelBasisSet", "alard-lupton")
policy.set("usePcaForSpatialKernel", True)
policy.set("spatialKernelOrder", 1)

spatialKernel, spatialBg, kernelCellSet = ipDiffim.createPsfMatchingKernel(imageToConvolve,
                                                                           imageToNotConvolve,
                                                                           policy)
 
cMi = afwImage.MaskedImageF(imageToConvolve.getDimensions())
afwMath.convolve(cMi, imageToConvolve, spatialKernel, False)
cMi.writeFits(outputImage)

