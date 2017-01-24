from __future__ import absolute_import, division, print_function

# all the c++ level classes and routines
from . import diffimLib

# all the other LSST packages
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

# Basically deprecated until SDQA is replaced


def makeRatingVector(kernelCellSet, spatialKernel, spatialBg):
    imstats = diffimLib.ImageStatisticsF()
    # sdqaVector = sdqa.SdqaRatingSet()

    width, height = spatialKernel.getDimensions()
    kImage = afwImage.ImageD(width, height)
    # find the kernel sum and its Rms by looking at the 4 corners of the image
    kSums = afwMath.vectorD()
    for x in (0, width):
        for y in (0, height):
            kSum = spatialKernel.computeImage(kImage, False, x, y)
            kSums.append(kSum)

    # afwStat    = afwMath.makeStatistics(kSums, afwMath.MEAN | afwMath.STDEV)
    # kSumRating = sdqa.SdqaRating("lsst.ip.diffim.kernel_sum",
    #                              afwStat.getValue(afwMath.MEAN),
    #                              afwStat.getValue(afwMath.STDEV),
    #                              scope)
    # sdqaVector.append(kSumRating)

    nGood = 0
    nBad = 0
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False):  # False = include bad candidates
            if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                # this has been used for processing
                nGood += 1

                xCand = int(cand.getXCenter())
                yCand = int(cand.getYCenter())

                # evaluate kernel and background at position of candidate
                kSum = spatialKernel.computeImage(kImage, False, xCand, yCand)
                kernel = afwMath.FixedKernel(kImage)
                background = spatialBg(xCand, yCand)

                diffIm = cand.getDifferenceImage(kernel, background)
                imstats.apply(diffIm)

                # candMean   = imstats.getMean()
                # candRms    = imstats.getRms()
                # candRating = sdqa.SdqaRating("lsst.ip.diffim.residuals_%d_%d" % (xCand, yCand),
                #                              candMean, candRms, scope)
                # sdqaVector.append(candRating)
            elif cand.getStatus() == afwMath.SpatialCellCandidate.BAD:
                nBad += 1

    # nGoodRating = sdqa.SdqaRating("lsst.ip.diffim.nCandGood", nGood, 0, scope)
    # sdqaVector.append(nGoodRating)
    # nBadRating = sdqa.SdqaRating("lsst.ip.diffim.nCandBad", nBad, 0, scope)
    # sdqaVector.append(nBadRating)

    nKernelTerms = spatialKernel.getNSpatialParameters()
    if nKernelTerms == 0:  # order 0
        nKernelTerms = 1
    # nBgTerms     = len(spatialBg.getParameters())
    # nKernRating  = sdqa.SdqaRating("lsst.ip.diffim.nTermsSpatialKernel", nKernelTerms, 0, scope)
    # nBgRating    = sdqa.SdqaRating("lsst.ip.diffim.nTermsSpatialBg", nBgTerms, 0, scope)
    # sdqaVector.append(nKernRating)
    # sdqaVector.append(nBgRating)

    # for i in range(sdqaVector.size()):
    #     pexLog.Trace("lsst.ip.diffim.makeSdqaRatingVector", 5,
    #                  "Sdqa Rating %s : %.2f %.2f" % (sdqaVector[i].getName(),
    #                                                  sdqaVector[i].getValue(),
    #                                                  sdqaVector[i].getErr()))
    #
    # return sdqaVector
