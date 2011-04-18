# all the c++ level classes and routines
import diffimLib

# all the other LSST packages
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math.mathLib as afwMath
import lsst.pex.logging as pexLog

# third party
import numpy
import time

#######
# Add noise
#######

def makeFlatNoiseImage(mi, seedStat = afwMath.MAX):
    img       = mi.getImage()
    seed      = int(10. * afwMath.makeStatistics(mi.getImage(), seedStat).getValue() + 1)
    rdm       = afwMath.Random(afwMath.Random.MT19937, seed)
    rdmImage  = img.Factory(img.getDimensions())
    afwMath.randomGaussianImage(rdmImage, rdm)
    return rdmImage

def makePoissonNoiseImage(im):
    """Return a Poisson noise image based on im
    
    Uses numpy.random; you may wish to call numpy.random.seed first.
    
    @warning This uses an undocumented numpy API (the documented API
    uses a single float expectation value instead of an array).
    
    @param[in] im image; the output image has the same dimensions and shape
        and its expectation value is the value of im at each pixel
    """
    import numpy.random as rand
    imArr = im.getArray()
    noiseIm = im.Factory(im.getBBox(afwImage.PARENT))
    noiseArr = noiseIm.getArray()

    intNoiseArr = rand.poisson(imArr)
    noiseArr[:, :] = intNoiseArr.astype(noiseArr.dtype)
    return noiseIm

#######
# Make fake images for testing; one is a delta function (or narrow
# gaussian) and the other is a convolution of this with a spatially
# varying kernel.
#######
def fakeCoeffs():
    kCoeffs = ((  6.9843,  0.0,         0.0), 
               (  3.2799, -0.0014,      0.0001), 
               (  4.6735,  0.0004,      6.3735e-06), 
               ( -6.4267, -0.0034,      0.0022), 
               ( -5.7412,  0.0015,      0.0009), 
               ( -4.0961, -0.0051,      0.0026), 
               ( -3.4660,  0.0003,      0.0007), 
               (  1.4134, -1.4700e-06, -0.0003), 
               (  0.9112, -0.0002,     -7.3208e-05), 
               ( -6.6402, -0.0025,      0.0003), 
               ( -5.4063, -0.0067,      0.0049), 
               (  2.7299,  0.0007,     -0.0008), 
               ( -3.5498, -0.0027,      0.0029), 
               ( -2.1196, -0.0013,      2.8696e-05), 
               ( -3.4298, -0.0065,      0.0053), 
               ( 37.9685,  0.0415,     -0.0314), 
               (  0.3578,  0.0009,     -0.0008), 
               (  1.7980,  0.0027,     -0.0002), 
               (-15.7181, -0.0176,      0.0129), 
               (  0.3994, -0.0005,      0.0003), 
               (-16.4162, -0.0174,      0.0125), 
               (  0.1032, -0.0002,      0.0002), 
               ( -0.3409, -0.0005,      6.1906e-05), 
               ( -0.3338, -0.0002,      0.0002), 
               ( -0.3492, -0.0011,      0.0001), 
               ( 14.6390,  0.0176,     -0.0119), 
               (  0.0118, -3.5455e-06, -1.2296e-05), 
               (  0.0662,  0.0002,     -8.8548e-05), 
               ( -2.9271, -0.0036,      0.0023), 
               (  0.0057,  0.0001,     -6.1779e-05), 
               ( -2.8832, -0.0033,      0.0022))
    return kCoeffs

def makeFakeKernelSet(policy, basisList, nCell = 5, deltaFunctionCounts = 1.e4, tGaussianWidth = 1.0,
                      addNoise = False, bgValue = 100., display = True):
    kSize    = policy.get('kernelSize')
    sizeCell = policy.get('sizeCellX')
    assert(sizeCell == policy.get('sizeCellY'))

    # This sets the final extent of each convolved delta function
    gaussKernelWidth   = sizeCell // 2

    # This sets the scale over which pixels are correlated in the
    # spatial convolution; should be at least as big as the kernel you
    # are trying to fit for
    spatialKernelWidth = kSize

    # Number of bad pixels due to convolutions
    border = (gaussKernelWidth + spatialKernelWidth)//2

    # Make a fake image with a matrix of delta functions
    totalSize = nCell * sizeCell + 2 * border
    tim       = afwImage.ImageF(totalSize, totalSize)
    for x in range(nCell):
        for y in range(nCell):
            tim.set(x * sizeCell + sizeCell // 2 + border - 1,
                    y * sizeCell + sizeCell // 2 + border - 1,
                    deltaFunctionCounts)

    # Turn this into stars with a narrow width; conserve counts
    gaussFunction = afwMath.GaussianFunction2D(tGaussianWidth, tGaussianWidth)
    gaussKernel   = afwMath.AnalyticKernel(gaussKernelWidth, gaussKernelWidth, gaussFunction)
    cim = afwImage.ImageF(tim.getDimensions())
    afwMath.convolve(cim, tim, gaussKernel, True)
    tim = cim

    # Trim off border pixels
    bbox = gaussKernel.shrinkBBox(tim.getBBox(afwImage.LOCAL))
    tim      = afwImage.ImageF(tim, bbox)
    # An estimate of its variance is itself
    tvar     = afwImage.ImageF(tim, True)
    # No mask
    tmask = afwImage.MaskU(tim.getDimensions())
    tmask.set(0x0)
    
    # Now make a science image which is this convolved with some
    # spatial function.  Use input basis list.
    #
    # THIS SHOULD BE SOMETHING WHOSE FIRST COMPONENT DOES NOT VARY
    # SPATIALLY, AND WHOSE OTHER COMPONENTS HAVE ZERO SUM.
    sOrder   = policy.get('spatialKernelOrder')
    polyFunc = afwMath.PolynomialFunction2D(sOrder)
    nParams  = len(polyFunc.getParameters())
    #kCoeffs  = []
    ## First one does not vary spatially
    #kCoeffs.append([])
    #kCoeffs[0].append(1.)
    #for i in range(1, nParams):
    #    kCoeffs[0].append(0.)
    #for i in range(1, len(basisList)):
    #    kCoeffs.append([])
    #    for j in range(nParams):
    #        kCoeffs[i].append(0.001 * (1.5 * (-1)**j + i))
    kCoeffs = fakeCoeffs()

    # Estimate of the image variance comes from convolution with main gaussian
    svar = afwImage.ImageF(tim.getDimensions())
    afwMath.convolve(svar, tim, basisList[0], False)
    # No mask
    smask = afwImage.MaskU(svar.getDimensions())
    smask.set(0x0)
    
    # Make the full convolved image
    sKernel = afwMath.LinearCombinationKernel(basisList, polyFunc)
    sKernel.setSpatialParameters(kCoeffs)
    cim = afwImage.ImageF(tim.getDimensions())
    afwMath.convolve(cim, tim, sKernel, False)

    # Get the good subregion
    bbox = sKernel.shrinkBBox(cim.getBBox(afwImage.LOCAL))

    # Add noise?
    if addNoise:
        svar += bgValue
        snoi  = makePoissonNoiseImage(svar)
        cim   = snoi
        cim  -= bgValue

        # noiseless?
        #tvar += bgValue
        #tnoi  = makePoissonNoiseImage(tvar)
        #tim  += tnoi

    # And turn into MaskedImages
    sim   = afwImage.ImageF(cim, bbox)
    svar  = afwImage.ImageF(svar, bbox)
    smask = afwImage.MaskU(smask, bbox)
    sMi   = afwImage.MaskedImageF(sim, smask, svar)
    
    tim   = afwImage.ImageF(tim, bbox)
    tvar  = afwImage.ImageF(tvar, bbox)
    tmask = afwImage.MaskU(tmask, bbox)
    tMi   = afwImage.MaskedImageF(tim, tmask, tvar)

    sMi.writeFits('science.fits')
    tMi.writeFits('template.fits')

    # look at this
    #import pylab
    #simarr = vectorFromImage(sim)
    #svararr = vectorFromImage(svar)
    #ssig = simarr / numpy.sqrt(svararr)
    #idx  = numpy.where(ssig < 10.)
    #ssig = ssig[idx]
    #pylab.figure()
    #n1, b1, p1 = pylab.hist(ssig, bins=100, normed=True)
    #ax2 = pylab.twinx()
    #ax2.plot(0.5 * (b1[1:] + b1[:-1]), n1)

    # no noise in template
    #timarr = vectorFromImage(tim)
    #tvararr = vectorFromImage(tvar)
    #tsig = timarr / numpy.sqrt(tvararr)
    #pylab.figure()
    #n2, b2, p2 = pylab.hist(tsig, bins=100, normed=True)
    #ax2 = pylab.twinx()
    #ax2.plot(0.5 * (b2[1:] + b2[:-1]), n2)

    #pylab.show()
    

    # No negative variance, which can screw things up
    assert(afwMath.makeStatistics(sMi.getVariance(), afwMath.MIN).getValue() >= 0.0)
    assert(afwMath.makeStatistics(tMi.getVariance(), afwMath.MIN).getValue() >= 0.0)

    # Get rid of any coordinate funniness
    sMi.setXY0(afwGeom.Point2I(0,0))
    tMi.setXY0(afwGeom.Point2I(0,0))

    if display:
        import lsst.afw.display.ds9 as ds9
        ds9.mtv(tMi.getImage(), frame=1)
        ds9.mtv(tMi.getVariance(), frame=2)
        ds9.mtv(sMi.getImage(), frame=3)
        ds9.mtv(sMi.getVariance(), frame=4)

    # Finally, make a kernelSet from these 2 images
    kernelCellSet = afwMath.SpatialCellSet(afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                                         afwGeom.Extent2I(sizeCell * nCell,
                                                                          sizeCell * nCell)),
                                           sizeCell,
                                           sizeCell)
    stampHalfWidth = 2 * kSize
    for x in range(nCell):
        for y in range(nCell):
            xCoord = x * sizeCell + sizeCell // 2
            yCoord = y * sizeCell + sizeCell // 2
            p0 = afwGeom.Point2I(xCoord - stampHalfWidth,
                                 yCoord - stampHalfWidth)
            p1 = afwGeom.Point2I(xCoord + stampHalfWidth,
                                 yCoord + stampHalfWidth)
            bbox = afwGeom.Box2I(p0, p1)
            tsi = afwImage.MaskedImageF(tMi, bbox, afwImage.LOCAL)
            ssi = afwImage.MaskedImageF(sMi, bbox, afwImage.LOCAL)

            kc = diffimLib.makeKernelCandidate(xCoord, yCoord, tsi, ssi, policy)
            kernelCellSet.insertCandidate(kc)

            #if display:
            #    ds9.mtv(tsi, frame=5)
            #    ds9.mtv(ssi, frame=6)
    
    return tMi, sMi, sKernel, kernelCellSet
    

#######
# Background subtraction for ip_diffim
#######

def backgroundSubtract(policy, maskedImages):
    t0 = time.time()
    algorithm   = policy.get("algorithm")
    binsize     = policy.get("binsize")
    undersample = policy.get("undersample")
    bctrl       = afwMath.BackgroundControl(algorithm)
    bctrl.setUndersampleStyle(undersample)
    for maskedImage in maskedImages:
        bctrl.setNxSample(maskedImage.getWidth()//binsize + 1)
        bctrl.setNySample(maskedImage.getHeight()//binsize + 1)
        
        image   = maskedImage.getImage() 
        backobj = afwMath.makeBackground(image, bctrl)
        image  -= backobj.getImageF()
        del image
        del backobj
        
    t1 = time.time()
    pexLog.Trace("lsst.ip.diffim.backgroundSubtract", 1,
                 "Total time for background subtraction : %.2f s" % (t1-t0))
    

#######
# Visualization of kernels
#######

def displayKernelList(kernelList, frame0 = 0):
    import lsst.afw.display.ds9 as ds9
    frame = frame0
    for kernel in kernelList:
        ki = afwImage.ImageD(kernel.getDimensions())
        kernel.computeImage(ki, False)
        ds9.mtv(ki, frame = frame)
        frame += 1

def displaySpatialKernelQuality(kernelCellSet, spatialKernel, spatialBg, frame):
    import lsst.afw.display.ds9 as ds9

    imstat = diffimLib.ImageStatisticsF()
    
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(True): # False = include bad candidates
            cand  = diffimLib.cast_KernelCandidateF(cand)
            if not (cand.getStatus() == afwMath.SpatialCellCandidate.GOOD):
                continue

            # raw inputs
            tmi = cand.getMiToConvolvePtr()
            smi = cand.getMiToNotConvolvePtr()
            ki  = cand.getImage()
            dmi = cand.getDifferenceImage(diffimLib.KernelCandidateF.RECENT)

            # spatial model
            ski   = afwImage.ImageD(ki.getDimensions())
            spatialKernel.computeImage(ski, False,
                                       afwImage.indexToPosition(int(cand.getXCenter())),
                                       afwImage.indexToPosition(int(cand.getYCenter())))
            sk    = afwMath.FixedKernel(ski)
            sbg   = spatialBg(afwImage.indexToPosition(int(cand.getXCenter())),
                              afwImage.indexToPosition(int(cand.getYCenter())))
            sdmi  = cand.getDifferenceImage(sk, sbg)

            ds9.mtv(tmi,  frame=frame+0) # template image
            ds9.mtv(smi,  frame=frame+1) # science image
            ds9.mtv(ki,   frame=frame+2) # best-fit kernel
            ds9.mtv(dmi,  frame=frame+3) # best-fit kernel diffim 
            ds9.mtv(ski,  frame=frame+4) # spatial kernel
            ds9.mtv(sdmi, frame=frame+5) # spatial kernel diffim

            ki -= ski
            ds9.mtv(ki,   frame=frame+6) # differnce in kernels

            imstat.apply(dmi)
            pexLog.Trace("lsst.ip.diffim.displaySpatialKernelQuality", 1,
                         "Candidate %d diffim residuals = %.2f +/- %.2f sigma" % (cand.getId(),
                                                                                  imstat.getMean(),
                                                                                  imstat.getRms()))
            imstat.apply(sdmi)
            pexLog.Trace("lsst.ip.diffim.displaySpatialKernelQuality", 1,
                         "Candidate %d sdiffim residuals = %.2f +/- %.2f sigma" % (cand.getId(),
                                                                                   imstat.getMean(),
                                                                                   imstat.getRms()))
            raw_input("Next: ")

                    
            
def displayKernelMosaic(kernelCellSet, frame):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.display.utils as displayUtils

    mos = displayUtils.Mosaic()
    
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False): # False = include bad candidates
            cand = diffimLib.cast_KernelCandidateF(cand)
            if not (cand.getStatus() == afwMath.SpatialCellCandidate.GOOD):
                continue
            im = cand.getKernelImage(diffimLib.KernelCandidateF.ORIG)
            mos.append(im)
            
    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)

def displaySpatialKernelMosaic(spatialKernel, width, height, frame, doNorm=False):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.display.utils as displayUtils

    mos = displayUtils.Mosaic()
    
    for y in (0, height//2, height):
        for x in (0, width//2, width):
            im   = afwImage.ImageD(spatialKernel.getDimensions())
            ksum = spatialKernel.computeImage(im,
                                              doNorm,
                                              afwImage.indexToPosition(x),
                                              afwImage.indexToPosition(y))
            mos.append(im, "x=%d y=%d kSum=%.2f" % (x, y, ksum))
    
    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)
    

def displayBasisMosaic(spatialKernel, frame):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.display.utils as displayUtils

    mos = displayUtils.Mosaic()

    basisList = spatialKernel.getKernelList()
    for idx in range(len(basisList)):
        kernel = basisList[idx]
        im   = afwImage.ImageD(spatialKernel.getDimensions())
        kernel.computeImage(im, False)
        mos.append(im, "K%d" % (idx))
    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)

def displayCandidateMosaic(kernelCellSet, frame):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.display.utils as displayUtils

    mos = displayUtils.Mosaic()

    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False): # False = include bad candidates
            cand  = diffimLib.cast_KernelCandidateF(cand)
            rchi2 = cand.getChi2()
                
            try:
                im = cand.getImage()
                if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    statStr = "Good"
                elif cand.getStatus() == afwMath.SpatialCellCandidate.BAD:
                    statStr = "Bad"
                else:
                    statStr = "Unkn"
                mos.append(im, "#%d: %.1f (%s)" % (cand.getId(), rchi2, statStr))
            except Exception, e:
                pass
    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)

def displayCandiateResults(kernelCellSet, frame):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.display.utils as displayUtils
    mos = displayUtils.Mosaic()

    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(True): # False = include bad candidates
            cand  = diffimLib.cast_KernelCandidateF(cand)
            try:
                tmi   = cand.getMiToConvolvePtr()
                smi   = cand.getMiToNotConvolvePtr()
                ki    = cand.getImage()
                dmi   = cand.returnDifferenceImage()

                mos.append(tmi)
                mos.append(foo1)
                #mos.append(ki)
                mos.append(foo2)
                
            except Exception, e:
                pass
            
    mosaic = mos.makeMosaic(mode=3)
    ds9.mtv(mosaic, frame=frame)

def displayFootprints(image, footprintList, frame):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.detection as afwDetection

    bitmask = image.getMask().getPlaneBitMask("DETECTED")
    afwDetection.setMaskFromFootprintList(image.getMask(), footprintList, bitmask)

    ds9.mtv(image, frame=frame)

    #ds9.mtv(image, frame=frame)
    #for fp in footprintList:
    #    bboxes = afwDetection.footprintToBBoxList(fp)
    #    for bbox in bboxes:
    #        bbox.shift(-afwGeom.Extent2I(image.getXY0()))
    #        
    #        x0, y0, x1, y1 = bbox.getMinX(), bbox.getMinY(), bbox.getMaxX(), bbox.getMaxY()
    #        
    #        x0 -= 0.5; y0 -= 0.5
    #        x1 += 0.5; y1 += 0.5
    #        
    #        ds9.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype=ds9.RED)
    

def displayKernelCellSet(image, kernelCellSet, frame):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.display.utils as utils
    
    ds9.ds9Cmd("mask transparency 50")
    ds9.mtv(image, frame = frame)

    for cell in kernelCellSet.getCellList():
        bbox = cell.getBBox()
        utils.drawBBox(bbox, frame = frame, ctype = 'blue')
        
        for cand in cell.begin(False): # False = include bad candidates
            cand  = diffimLib.cast_KernelCandidateF(cand)
   
            if (cand.getStatus() == afwMath.SpatialCellCandidate.GOOD):
                color = 'green'
            elif (cand.getStatus() == afwMath.SpatialCellCandidate.BAD):
                color = 'red'
            else:
                color = 'yellow'
            
            xCenter = cand.getXCenter()
            yCenter = cand.getYCenter()

            cmd   = 'regions command { circle %g %g %g # color=%s text=\"%d\"};' % (
                xCenter, yCenter, 10, color, cand.getId())
            ds9.ds9Cmd(cmd)
