# all the c++ level classes and routines
import diffimLib

# other diffim routines
from rejectKernelSumOutliers import rejectKernelSumOutliers
from runPca import runPca
from diffimTools import findIqr, vectorFromImage, imageFromVector

# all the other LSST packages
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9

# third party
import numpy as num

display = False

def returnMeanKernel(templateMI, scienceMI, footprintList, policy, growFp=False, useAlard=False, renormalize=True):
    # template image is to be convolved
    # footprint list is not assumed to be clean, and to have come straight from detection
    # if you have *already* grown the Footprints, set growFp=False
    # if you want to enfore all kernels to have the same (mean) kernel sum use renormalize=True

    # Same dimensions, non negotiable
    assert (templateMI.getDimensions() == scienceMI.getDimensions())

    # We also assume that at this stage, they are aligned at the pixel level
    # Assign to the coordinate system of the science image
    templateMI.setXY0( scienceMI.getXY0() )

    # Create PSF matching functor
    kCols = policy.get("kernelCols")
    kRows = policy.get("kernelRows")
    if useAlard:
        nGauss   = policy.get("alardNGauss")
        sigGauss = policy.getDoubleArray("alardSigGauss")
        degGauss = policy.getIntArray("alardDegGauss")
        assert len(sigGauss) == nGauss
        assert len(degGauss) == nGauss
        assert kCols == kRows  # square
        assert kCols % 2 == 1  # odd sized
        kHalfWidth = int(kCols/2)
        kBasisList = diffimLib.generateAlardLuptonKernelSet(kHalfWidth, nGauss, sigGauss, degGauss)
        kFunctor   = diffimLib.PsfMatchingFunctorF(kBasisList)
    else:
        kBasisList = diffimLib.generateDeltaFunctionKernelSet(kCols, kRows)
        if policy.get("regularizationUse") == True:
            order      = policy.get("regularizationOrder")
            boundary   = policy.get("regularizationBoundary")
            difference = policy.get("regularizationDifference")
            H          = diffimLib.generateFiniteDifferenceRegularization(kCols, kRows, order, boundary, difference)
            kFunctor   = diffimLib.PsfMatchingFunctorF(kBasisList, H)
        else:
            kFunctor = diffimLib.PsfMatchingFunctorF(kBasisList)

    # Look for set bits
    itcFunctor  = diffimLib.FindSetBitsU(templateMI.getMask()) 
    itncFunctor = diffimLib.FindSetBitsU(scienceMI.getMask())

    # Kernel quality
    dStats          = diffimLib.ImageStatisticsF()
    maxResidualMean = policy.get("maximumFootprintResidualMean")
    maxResidualStd  = policy.get("maximumFootprintResidualStd")

    # How much you grow each raw detection
    fpGrowKsize = policy.get("fpGrowKsize")
    fpGrowPix   = int(fpGrowKsize * kCols)

    # List of good kernels
    kernelList = []
    for fp in footprintList:
        if growFp:
            fpUse = afwDetection.growFootprint(fp, fpGrowPix, False)
        else:
            fpUse = fp
            
        bbox   = fpUse.getBBox()
        bbox.shift(-scienceMI.getX0(), -scienceMI.getY0());

        # grown off the image?
        if ( (bbox.getX0() < 0) or
             (bbox.getY0() < 0) or
             (bbox.getX1() > scienceMI.getWidth()) or
             (bbox.getY1() > scienceMI.getHeight()) ):
            continue

        # check for bad pixels
        itcFunctor.apply(fpUse)
        itncFunctor.apply(fpUse)

        if ( itcFunctor.getBits() or itncFunctor.getBits() ):
            continue

        # looks like a good kernel
        tmi  = afwImage.MaskedImageF(templateMI, bbox)
        smi  = afwImage.MaskedImageF(scienceMI, bbox)
        # variance estimate
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        if display:
            ds9.mtv(tmi, frame=1)
            ds9.mtv(smi, frame=0)

        # solve for kernel
        kFunctor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), policy)
        kernel    = kFunctor.getKernel()
        kImageOut = afwImage.ImageD(kCols, kRows)
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = diffimLib.convolveAndSubtract(tmi, smi, kernel, kFunctor.getBackground())
        dStats.apply( diffIm )

        if display:
            ds9.mtv(kImageOut, frame=2)
            ds9.mtv(diffIm,    frame=3)

        print 'DF Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; ' % (dStats.getMean(), dStats.getRms(),
                                                                     kSum, kFunctor.getBackground()),
        
        if (abs(dStats.getMean()) > maxResidualMean) or (dStats.getRms() > maxResidualStd):
            print 'Rejected'
            continue
        print 'OK'
        kernelList.append( (kSum, kernel, kImageOut) )

    # Do rejection of kSum outliers
    maxOutlierIterations = policy.get('maxOutlierIterations')
    maxOutlierSigma      = policy.get('maxOutlierSigma')
    newList              = kernelList
    for nIter in range(maxOutlierIterations):
        kSums  = num.array( [x[0] for x in newList] )
        kSumMedian, kSumIqr, kSumSigMean, kSumSigMedian = findIqr(kSums)
        print 'Niter %d, N = %d : %.3f %.3f %.3f %.3f' % (nIter, len(newList), kSumMedian, kSumIqr, kSumSigMean, kSumSigMedian)
        goodIdx = num.where ( (kSums-kSumMedian) < (kSumSigMedian*maxOutlierSigma) )[0]
        if len(goodIdx) == len(newList):
            break
        else:
            kernelList = []
            for idx in goodIdx:
                kernelList.append( newList[idx] )
            newList = kernelList
    # Upon exit
    kernelList = newList
    nKernels = len(kernelList)

    # Renormalize
    if renormalize:
        kernelListNew = []
        for i in range(nKernels):
            kSum, kernel, kImageOut = kernelList[i]
            kImageOut              /= kSum        # give kernel sum of 1
            kImageOut              *= kSumMedian  # give kernel sum of kSumMedian
            kernel                  = afwMath.FixedKernel(kImageOut)
            kernelListNew.append( (kSumMedian, kernel, kImageOut) )
        kernelList = kernelListNew
        
    # Find basic mean shape
    M        = num.zeros((kCols*kRows, nKernels))
    for idx in range(nKernels):
        M[:,idx] = vectorFromImage(kernelList[idx][2])
            
    mean, U, eVal, eCoeff = runPca(M, None)
    meanImage  = imageFromVector(mean, kCols, kRows, retType=afwImage.ImageD)
    meanKernel = afwMath.FixedKernel(meanImage)

    if display:
        frame = 0
        for i in range(nKernels):
            ds9.mtv( kernelList[i][2], frame=frame )
            frame += 1
        ds9.mtv(meanImage, frame=frame)

    return meanKernel
