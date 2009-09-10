import numpy
import time

# all the c++ level classes and routines
import diffimLib

# all the other diffim routines
from createSpatialModelKernelCells import createSpatialModelKernelCells
from rejectKernelSumOutliers import rejectKernelSumOutliers
from spatialKernelFit import spatialModelKernelPca, spatialModelByPca, evaluateModelByPca, \
                             spatialModelByPixel, evaluateModelByPixel

# all the other LSST packages
import lsst.pex.logging as pexLog
import lsst.afw.display.ds9 as ds9
import lsst.afw.math  as afwMath
import lsst.afw.image as afwImage
import lsst.sdqa as sdqa

def subtractMaskedImage(templateMaskedImage, 
                        scienceMaskedImage, 
                        policy, 
                        log,
                        fpList=None,
                        useAlard=False,
                        invert=False,
                        display=False):
    
    # Make sure they are the EXACT same dimensions in pixels
    # This is non-negotiable
    assert (templateMaskedImage.getDimensions() == \
            scienceMaskedImage.getDimensions())
    
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
            order    = policy.get("regularizationOrder")
            H        = diffimLib.generateDeltaFunctionRegularization(kCols, kRows, order)
            kFunctor = diffimLib.PsfMatchingFunctorF(kBasisList, H)
        else:
            kFunctor = diffimLib.PsfMatchingFunctorF(kBasisList)
        
        

    if fpList == None:
        # Need to find own footprints
        log.log(pexLog.Log.INFO, "Starting footprints : %s" % (time.ctime()))

        fpList = diffimLib.getCollectionOfFootprintsForPsfMatching(
            templateMaskedImage,
            scienceMaskedImage,
            policy)
        log.log(pexLog.Log.INFO, "Ending footprints : %s" % (time.ctime()))

    if display:
        frame=1

        diffimPlaneName = "DIFFIM_STAMP_PLANE"
        diffimStampPlane = scienceMaskedImage.getMask().addMaskPlane(diffimPlaneName)
        ds9.setMaskPlaneColor(diffimPlaneName, "cyan") # afw > 3.3.10 can handle any X11 colour, e.g. "powderBlue"

        afwDetection.setMaskFromFootprintList(scienceMaskedImage.getMask(), fpList, (1 << diffimStampPlane))
        ds9.mtv(scienceMaskedImage.getMask(), frame=frame)
        
        for fp in fpList:
            x0, y0 = fp.getBBox().getLLC() - scienceMaskedImage.getXY0()
            x1, y1 = fp.getBBox().getURC() - scienceMaskedImage.getXY0()
            ds9.line([(x0, y0), (x0, y1), (x1, y1), (x1, y0), (x0, y0)], frame=frame)

    # switch image you convolve
    if invert:
        (scienceMaskedImage, templateMaskedImage) = (templateMaskedImage, scienceMaskedImage)

    # Set up grid for spatial model
    log.log(pexLog.Log.INFO, "Starting kernel : %s" % (time.ctime()))
    spatialCells = createSpatialModelKernelCells(
        templateMaskedImage,
        scienceMaskedImage,
        fpList,
        kFunctor,
        policy,
        display=display)

    # Reject kernel sum outliers
    rejectKernelSumOutliers(spatialCells, policy)

    # Set up fitting loop 
    maxSpatialIterations = policy.getInt("maxSpatialIterations")
    rejectKernels = policy.getBool("spatialKernelRejection")
    nRejected = -1
    nIter =  0
    
    # And fit spatial kernel model
    if policy.get("spatialKernelModel") == "pca":
        # Fit spatial variation of principal components

        minPrincipalComponents = policy.getInt("minPrincipalComponents")
        maxPrincipalComponents = policy.getInt("maxPrincipalComponents")
        fracEigenVal = policy.getDouble("fracEigenVal")
        
        while (nRejected != 0) and (nIter < maxSpatialIterations):
            # Run the PCA
            mKernel, eKernelVector, eVal, eCoeff = \
                    spatialModelKernelPca(spatialCells, policy)

            # Make the decision on how many components to use
            eFrac = numpy.cumsum(eVal)
            eFrac /= eFrac[-1]
            nEval = len(numpy.where(eFrac < fracEigenVal)[0])
            nEval = min(nEval, maxPrincipalComponents)
            nEval = max(nEval, minPrincipalComponents)

            pexLog.Trace("lsst.ip.diffim.subtractMaskedImage", 3, 
                         "PCA iteration %d : Using %d principal components" % (nIter, nEval))

            # do spatial fit here by Principal Component
            sKernel, bgFunction = spatialModelByPca(
                spatialCells,
                mKernel,
                eKernelVector,
                eCoeff,
                nEval,
                policy)

            # Evaluate quality of spatial fit
            nRejected, sdqaList = evaluateModelByPca(
                spatialCells,
                bgFunction, 
                sKernel,
                policy, 
                reject=rejectKernels)

            # If you get a new kernel, make sure its consistent with the kernel sums
            if nRejected:
                rejectKernelSumOutliers(spatialCells, policy)
                
            nIter += 1

    elif policy.get("spatialKernelModel") == "pixel":
        # Fit function to each pixel

        while (nRejected != 0) and (nIter < maxSpatialIterations):
            # do spatial fit here pixel by pixel
            sKernel, bgFunction = spatialModelByPixel(
                spatialCells,
                kBasisList,
                policy)

            # and check quality
            nRejected, sdqaList  = evaluateModelByPixel(
                spatialCells,
                bgFunction, 
                sKernel, 
                policy, 
                reject=rejectKernels)

            # If you get a new kernel, make sure its consistent with the kernel sums
            if nRejected:
                rejectKernelSumOutliers(spatialCells, policy)

            nIter += 1
        
    else:
        # All that is supported
        # Throw exception!
        pass

    log.log(pexLog.Log.INFO, "Ending kernel : %s" % (time.ctime()))
    log.log(pexLog.Log.INFO, "Starting convolve : %s" % (time.ctime()))
    if policy.exists("backgroundPolicy"):
        background = 0                  # no need to subtract a background in subtraction as we'll do so in a moment
    else:
        background = bgFunction
        
    differenceMaskedImage = diffimLib.convolveAndSubtract(templateMaskedImage,
                                                         scienceMaskedImage, sKernel, background)
    log.log(pexLog.Log.INFO, "Ending convolve : %s" % (time.ctime()))

    #
    # Maybe subtract a background model from the difference image
    #
    if policy.exists("backgroundPolicy"):
        algorithm = policy.get("backgroundPolicy.algorithm")
        binsize   = policy.get("backgroundPolicy.binsize")

        if algorithm == "NATURAL_SPLINE":
            bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
        else:
            raise RuntimeError, "Unknown backgroundPolicy.algorithm: %s" % (algorithm)

        bctrl.setNxSample(int(differenceMaskedImage.getWidth()//binsize) + 1)
        bctrl.setNySample(int(differenceMaskedImage.getHeight()//binsize) + 1)
        
        image = differenceMaskedImage.getImage() 
        backobj = afwMath.makeBackground(image, bctrl)
        image -= backobj.getImageF()
        del image; del backobj

    if display:
        frame = 3
        ds9.mtv(differenceMaskedImage, frame=frame)
        ds9.dot("Subtracted", 0, 0, frame=frame)

        if False:
            chisqMI = differenceMaskedImage.Factory(differenceMaskedImage, True)
            chisq = chisqMI.getImage();
            chisq *= chisq; chisq /= scienceMaskedImage.getVariance()
            del chisq

            frame = 4
            ds9.mtv(chisqMI, frame=frame)

    # If we inverted the difference image...
    if invert:
        differenceMaskedImage *= -1

    # N.b. Per-footprint sdqa ratings are not implemented for DC3a.
    # Override the list returned from evaluateModelBy... for now.
    sdqaList = sdqa.SdqaRatingSet()
    
    #
    # Lets do some more Sdqa here
    #
    imStats = diffimLib.ImageStatisticsF()
    imStats.apply(differenceMaskedImage)
    residualsRating = sdqa.SdqaRating("ip.diffim.residuals",
            imStats.getMean(), 
            imStats.getRms(), 
            sdqa.SdqaRating.AMP) 
    sdqaList.append(residualsRating)

    # Check kernel sum in the corners
    kSums  = []
    kImage = afwImage.ImageD(sKernel.getDimensions())
    for nRow in [0, templateMaskedImage.getHeight()]:
        for nCol in [0, templateMaskedImage.getWidth()]:
            kSums.append( sKernel.computeImage(kImage, False, nCol, nRow) )
    kSumArray = numpy.array(kSums)
    pexLog.Trace("lsst.ip.diffim.subtractMaskedImage", 3, 
            "Final Kernel Sum from Image Corners : %0.3f (%0.3f)" % 
            (kSumArray.mean(), kSumArray.std()))
    kernelSumRating =  sdqa.SdqaRating("ip.diffim.kernelSum",
            kSumArray.mean(),
            kSumArray.std(), 
            sdqa.SdqaRating.AMP)  
    sdqaList.append(kernelSumRating)

    # What kind of metadata do we add here to MaskedImage?
    
    return differenceMaskedImage, sKernel, bgFunction, sdqaList
