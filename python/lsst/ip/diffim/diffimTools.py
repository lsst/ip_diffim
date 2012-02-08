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

# all the other LSST packages
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math.mathLib as afwMath
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
from psfMatch import PsfMatchConfigAL
from makeKernelBasisList import makeKernelBasisList 

# third party
import numpy
import time
import os

# Helper functions for ipDiffim; mostly viewing of results and writing
# debugging info to disk.

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
    kCoeffs = ((  1.0,     0.0,       0.0), 
               (  0.005,  -0.000001,  0.000001), 
               (  0.005,   0.000004,  0.000004), 
               ( -0.001,  -0.000030,  0.000030), 
               ( -0.001,   0.000015,  0.000015), 
               ( -0.005,  -0.000050,  0.000050))
    return kCoeffs

def makeFakeKernelSet(sizeCell = 128, nCell = 3,
                      deltaFunctionCounts = 1.e4, tGaussianWidth = 1.0,
                      addNoise = True, bgValue = 100., display = False):

    configFake               = PsfMatchConfigAL()
    configFake.alardNGauss   = 1
    configFake.alardSigGauss = [2.5,]
    configFake.alardDegGauss = [2,]
    configFake.sizeCellX     = sizeCell
    configFake.sizeCellY     = sizeCell
    configFake.spatialKernelOrder = 1.0
    configFake.spatialModelType = "polynomial"
    configFake.singleKernelClipping = False   # variance is a hack
    configFake.spatialKernelClipping = False  # variance is a hack
    if bgValue > 0.0:
        configFake.fitForBackground = True

    policyFake = pexConfig.makePolicy(configFake)

    basisList = makeKernelBasisList(configFake)
    kSize     = configFake.kernelSize

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
    tim       = afwImage.ImageF(afwGeom.Extent2I(totalSize, totalSize))
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
    tim  = afwImage.ImageF(tim, bbox, afwImage.LOCAL)
    
    # Now make a science image which is this convolved with some
    # spatial function.  Use input basis list.
    sOrder   = 1
    polyFunc = afwMath.PolynomialFunction2D(1)
    kCoeffs  = fakeCoeffs()
    nToUse   = min(len(kCoeffs), len(basisList))

    # Make the full convolved science image
    sKernel = afwMath.LinearCombinationKernel(afwMath.KernelList(basisList[:nToUse]), polyFunc)
    sKernel.setSpatialParameters(kCoeffs[:nToUse])
    sim = afwImage.ImageF(tim.getDimensions())
    afwMath.convolve(sim, tim, sKernel, True)

    # Get the good subregion
    bbox = sKernel.shrinkBBox(sim.getBBox(afwImage.LOCAL))

    # Add background
    sim  += bgValue

    # Watch out for negative values
    tim  += 2 * numpy.abs(numpy.min(tim.getArray()))
    
    # Add noise?
    if addNoise:
        sim   = makePoissonNoiseImage(sim)
        tim   = makePoissonNoiseImage(tim) 
    
    # And turn into MaskedImages
    sim   = afwImage.ImageF(sim, bbox, afwImage.LOCAL)
    svar  = afwImage.ImageF(sim, True)
    svar.set(1.0)
    smask = afwImage.MaskU(sim.getDimensions())
    smask.set(0x0)
    sMi   = afwImage.MaskedImageF(sim, smask, svar)
    
    tim   = afwImage.ImageF(tim, bbox, afwImage.LOCAL)
    tvar  = afwImage.ImageF(tim, True)
    tvar.set(1.0)
    tmask = afwImage.MaskU(tim.getDimensions())
    tmask.set(0x0)
    tMi   = afwImage.MaskedImageF(tim, tmask, tvar)

    if display:
        import lsst.afw.display.ds9 as ds9
        ds9.mtv(tMi, frame=1)
        ds9.mtv(sMi, frame=2)

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

            kc = diffimLib.makeKernelCandidate(xCoord, yCoord, tsi, ssi, policyFake)
            kernelCellSet.insertCandidate(kc)

    return tMi, sMi, sKernel, kernelCellSet, configFake
    

#######
# Background subtraction for ip_diffim
#######

def backgroundSubtract(config, maskedImages):
    backgrounds = []
    t0 = time.time()
    algorithm   = config.algorithmName
    binsize     = config.binsize
    undersample = config.undersample
    bctrl       = afwMath.BackgroundControl(algorithm)
    bctrl.setUndersampleStyle(undersample)
    for maskedImage in maskedImages:
        bctrl.setNxSample(maskedImage.getWidth()//binsize + 1)
        bctrl.setNySample(maskedImage.getHeight()//binsize + 1)
        
        image   = maskedImage.getImage() 
        backobj = afwMath.makeBackground(image, bctrl)
        image  -= backobj.getImageF()
        backgrounds.append(backobj.getImageF())
        del image
        del backobj
        
    t1 = time.time()
    pexLog.Trace("lsst.ip.diffim.backgroundSubtract", 1,
                 "Total time for background subtraction : %.2f s" % (t1-t0))
    return backgrounds

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
                                       int(cand.getXCenter()),
                                       int(cand.getYCenter()))
            sk    = afwMath.FixedKernel(ski)
            sbg   = spatialBg(int(cand.getXCenter()),
                              int(cand.getYCenter()))
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
            ksum = spatialKernel.computeImage(im, doNorm, x, y)

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

def displayCandidateResults(kernelCellSet, frame, goodOnly = True):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.display.utils as displayUtils
    mos = displayUtils.Mosaic()

    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(goodOnly): # False = include bad candidates
            cand  = diffimLib.cast_KernelCandidateF(cand)
            
            tmi   = cand.getMiToConvolvePtr()
            smi   = cand.getMiToNotConvolvePtr()
            try:
                ki    = cand.getKernelImage(diffimLib.KernelCandidateF.ORIG)
                dmi   = cand.getDifferenceImage(diffimLib.KernelCandidateF.ORIG)
            except:
                pass
            else:
                mos.append(tmi.getImage())
                mos.append(smi.getImage())
                mos.append(ki.convertFloat())
                mos.append(dmi.getImage())
                
    mosaic = mos.makeMosaic(mode=4)
    ds9.mtv(mosaic, frame=frame)

def displayFootprints(image, footprintList, frame):
    import lsst.afw.display.ds9 as ds9
    import lsst.afw.detection as afwDetection

    bitmask = image.getMask().getPlaneBitMask("DETECTED")
    afwDetection.setMaskFromFootprintList(image.getMask(), footprintList, bitmask)

    ds9.mtv(image, frame=frame)

    # ANOTHER WAY TO DO IT
    #
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
    

def writeKernelCellSet(kernelCellSet, psfMatchingKernel, backgroundModel, outdir):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False): # False = include bad candidates
            cand = diffimLib.cast_KernelCandidateF(cand)
            if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                xCand = int(cand.getXCenter())
                yCand = int(cand.getYCenter())
                idCand = cand.getId()
                diffIm = cand.getDifferenceImage(diffimLib.KernelCandidateF.ORIG)
                kernel = cand.getKernelImage(diffimLib.KernelCandidateF.ORIG)
                diffIm.writeFits(os.path.join(outdir, 'diffim_c%d_x%d_y%d.fits' % (idCand, xCand, yCand)))
                kernel.writeFits(os.path.join(outdir, 'kernel_c%d_x%d_y%d.fits' % (idCand, xCand, yCand)))

                # Diffim from spatial model
                ski   = afwImage.ImageD(kernel.getDimensions())
                psfMatchingKernel.computeImage(ski, False, xCand, yCand)
                sk    = afwMath.FixedKernel(ski)
                sbg   = backgroundModel(xCand, yCand)
                sdmi  = cand.getDifferenceImage(sk, sbg)
                sdmi.writeFits(os.path.join(outdir, 'sdiffim_c%d_x%d_y%d.fits' % (idCand, xCand, yCand)))
                                
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
