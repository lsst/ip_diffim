# all the c++ level classes and routines
import diffimLib

# all the other LSST packages
import lsst.afw.image as afwImage
import lsst.afw.math.mathLib as afwMath
import lsst.pex.logging as pexLog

# third party
import numpy
import time

#######
# Expansions of functionality found in lsst.afw.image.testUtils
#######

def vectorFromImage(im, dtype=float):
    vec = numpy.zeros(im.getWidth()*im.getHeight(), dtype=dtype)
    idx = 0
    for row in range(im.getHeight()):
        for col in range(im.getWidth()):
            vec[idx] = im.get(col, row)
            idx     += 1
    return vec

def imageFromVector(vec, width, height, retType=afwImage.ImageF):
    im  = retType(width, height)
    idx = 0
    for row in range(height):
        for col in range(width):
            # need to cast numpy.float64 as float
            im.set(col, row, float(vec[idx]))
            idx     += 1
    return im

#######
# Background subtraction for ip_diffim
#######

def backgroundSubtract(policy, maskedImages):
    t0 = time.time()
    algorithm   = policy.get("backgroundPolicy.algorithm")
    binsize     = policy.get("backgroundPolicy.binsize")
    undersample = policy.get("backgroundPolicy.undersample")
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
            dmi = cand.returnDifferenceImage()

            # spatial model
            ski   = afwImage.ImageD(ki.getDimensions())
            spatialKernel.computeImage(ski, False,
                                       afwImage.indexToPosition(int(cand.getXCenter())),
                                       afwImage.indexToPosition(int(cand.getYCenter())))
            sk    = afwMath.FixedKernel(ski)
            sbg   = spatialBg(afwImage.indexToPosition(int(cand.getXCenter())),
                              afwImage.indexToPosition(int(cand.getYCenter())))
            sdmi  = cand.returnDifferenceImage(sk, sbg)

            ds9.mtv(tmi,  frame=frame+0)
            ds9.mtv(smi,  frame=frame+1)
            ds9.mtv(ki,   frame=frame+2)
            ds9.mtv(dmi,  frame=frame+3)
            ds9.mtv(ski,  frame=frame+4)
            ds9.mtv(sdmi, frame=frame+5)

            ki -= ski
            ds9.mtv(ki,   frame=frame+6)

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
            #raw_input("Next: ")

                    
            
def displaySpatialKernelMosaic(spatialKernel, width, height, frame):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.display.utils as displayUtils

    mos = displayUtils.Mosaic()
    
    for x in (0, width//2, width):
        for y in (0, height//2, height):
            im   = afwImage.ImageD(spatialKernel.getDimensions())
            ksum = spatialKernel.computeImage(im,
                                              False,
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
                mos.append(smi)
                mos.append(ki)
                mos.append(dmi)
                
            except Exception, e:
                pass
            
    mosaic = mos.makeMosaic(mode=4)
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
    #        bbox.shift(-image.getX0(),-image.getY0())
    #        
    #        x0, y0, x1, y1 = bbox.getX0(), bbox.getY0(), bbox.getX1(), bbox.getY1()
    #        
    #        x0 -= 0.5; y0 -= 0.5
    #        x1 += 0.5; y1 += 0.5
    #        
    #        ds9.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype=ds9.RED)
    

