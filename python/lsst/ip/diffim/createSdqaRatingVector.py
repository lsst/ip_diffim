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

# all the c++ level classes and routines
import diffimLib

# all the other diffim routines
import lsst.sdqa as sdqa

# all the other LSST packages
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.logging as pexLog

def createSdqaRatingVector(kernelCellSet, spatialKernel, spatialBg, scope=sdqa.SdqaRating.AMP):
    imstats    = diffimLib.ImageStatisticsF()
    sdqaVector = sdqa.SdqaRatingSet()

    width, height = spatialKernel.getDimensions()
    kImage        = afwImage.ImageD(width, height)
    # find the kernel sum and its Rms by looking at the 4 corners of the image
    kSums = afwMath.vectorD()
    for x in (0, width):
        for y in (0, height):
            kSum = spatialKernel.computeImage(kImage, False,
                                              afwImage.indexToPosition(x),
                                              afwImage.indexToPosition(y))
            kSums.push_back(kSum)
            
    afwStat    = afwMath.makeStatistics(kSums, afwMath.MEAN | afwMath.STDEV)
    kSumRating = sdqa.SdqaRating("lsst.ip.diffim.kernel_sum",
                                 afwStat.getValue(afwMath.MEAN),
                                 afwStat.getValue(afwMath.STDEV),
                                 scope)
    sdqaVector.append(kSumRating)
    
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(True): # only look at non-bad candidates
            cand = diffimLib.cast_KernelCandidateF(cand)
            if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                # this has been used for processing
    
                xCand = int(cand.getXCenter())
                yCand = int(cand.getYCenter())

                # evaluate kernel and background at position of candidate
                kSum = spatialKernel.computeImage(kImage, False,
                                                  afwImage.indexToPosition(xCand),
                                                  afwImage.indexToPosition(yCand))
                kernel = afwMath.FixedKernel(kImage)
                
                background = spatialBg(afwImage.indexToPosition(xCand),
                                       afwImage.indexToPosition(yCand))

                diffIm = cand.returnDifferenceImage(kernel, background)
                imstats.apply(diffIm)
                
                candMean   = imstats.getMean()
                candRms    = imstats.getRms()
                candRating = sdqa.SdqaRating("lsst.ip.diffim.residuals_%d_%d" % (xCand, yCand),
                                             candMean, candRms, scope)
                sdqaVector.append(candRating)

    for i in range(sdqaVector.size()):
        pexLog.Trace("lsst.ip.diffim.createSdqaRatingVector", 3,
                     "Sdqa Rating %s : %.2f %.2f" % (sdqaVector[i].getName(),
                                                     sdqaVector[i].getValue(),
                                                     sdqaVector[i].getErr()))
                     
    return sdqaVector
