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

# python
import time
import os
from collections import Counter
import numpy as np

# all the c++ level classes and routines
import diffimLib

# all the other LSST packages
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDetect
import lsst.afw.math.mathLib as afwMath
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.meas.algorithms as measAlg
import lsst.meas.deblender.baseline as measDeblend
from .makeKernelBasisList import makeKernelBasisList 

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

    from . import imagePsfMatch
    configFake               = imagePsfMatch.ImagePsfMatchConfig()
    configFake.kernel.name   = "AL"
    subconfigFake            = configFake.kernel.active
    subconfigFake.alardNGauss   = 1
    subconfigFake.alardSigGauss = [2.5,]
    subconfigFake.alardDegGauss = [2,]
    subconfigFake.sizeCellX     = sizeCell
    subconfigFake.sizeCellY     = sizeCell
    subconfigFake.spatialKernelOrder = 1
    subconfigFake.spatialModelType = "polynomial"
    subconfigFake.singleKernelClipping = False   # variance is a hack
    subconfigFake.spatialKernelClipping = False  # variance is a hack
    if bgValue > 0.0:
        subconfigFake.fitForBackground = True

    policyFake = pexConfig.makePolicy(subconfigFake)

    basisList = makeKernelBasisList(subconfigFake)
    kSize     = subconfigFake.kernelSize

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
    tim  += 2 * np.abs(np.min(tim.getArray()))
    
    # Add noise?
    if addNoise:
        sim   = makePoissonNoiseImage(sim)
        tim   = makePoissonNoiseImage(tim) 
    
    # And turn into MaskedImages
    sim   = afwImage.ImageF(sim, bbox, afwImage.LOCAL)
    svar  = afwImage.ImageF(sim, True)
    smask = afwImage.MaskU(sim.getDimensions())
    smask.set(0x0)
    sMi   = afwImage.MaskedImageF(sim, smask, svar)
    
    tim   = afwImage.ImageF(tim, bbox, afwImage.LOCAL)
    tvar  = afwImage.ImageF(tim, True)
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
        del backobj

    t1 = time.time()
    pexLog.Trace("lsst.ip.diffim.backgroundSubtract", 1,
                 "Total time for background subtraction : %.2f s" % (t1-t0))
    return backgrounds

#######
# More coarse debugging
#######

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
                                
#######
# Converting types
#######

def sourceToFootprintList(candidateInList, templateExposure, scienceExposure, config, log):
    """ Takes an input list of Sources that were selected to constrain
    the Psf-matching Kernel and turns them into a List of Footprints,
    which are used to seed a set of KernelCandidates.  The function
    checks both the template and science image for masked pixels,
    rejecting the Source if certain Mask bits (defined in config) are
    set within the Footprint.

    @param candidateInList: Input list of Sources
    @param templateExposure: Template image, to be checked for Mask bits in Source Footprint
    @param scienceExposure: Science image, to be checked for Mask bits in Source Footprint
    @param config: Config that defines the Mask planes that indicate an invalid Source and Bbox grow radius
    @param log: Log for output
    
    @return a List of Footprints whose pixels will be used to constrain the Psf-matching Kernel
    """

    candidateOutList = []
    fsb = diffimLib.FindSetBitsU()
    badBitMask = 0
    for mp in config.badMaskPlanes: 
        badBitMask |= afwImage.MaskU.getPlaneBitMask(mp)
    log.info("Growing %d kernel candidate stars" % (len(candidateInList)))
    bbox = scienceExposure.getBBox()
    for kernelCandidate in candidateInList:
        if not type(kernelCandidate) == afwTable.SourceRecord:
            raise RuntimeError, ("Candiate not of type afwTable.SourceRecord")
        bm1 = 0
        bm2 = 0
        ra  = kernelCandidate.getRa()
        dec = kernelCandidate.getDec()
        center = afwGeom.Point2I(scienceExposure.getWcs().skyToPixel(kernelCandidate.getCoord()))
        if center[0] < bbox.getMinX() or center[0] > bbox.getMaxX():
            continue
        if center[1] < bbox.getMinY() or center[1] > bbox.getMaxY():
            continue
        # Grow Sources
        xmin   = center[0] - config.fpGrowPix
        xmax   = center[0] + config.fpGrowPix
        ymin   = center[1] - config.fpGrowPix
        ymax   = center[1] + config.fpGrowPix

        # Keep object centered
        if (xmin - bbox.getMinX()) < 0:
            xmax += (xmin - bbox.getMinX())
            xmin -= (xmin - bbox.getMinX())
        if (ymin - bbox.getMinY()) < 0:
            ymax += (ymin - bbox.getMinY())
            ymin -= (ymin - bbox.getMinY())
        if (bbox.getMaxX() - xmax) < 0:
            xmin -= (bbox.getMaxX() - xmax)
            xmax += (bbox.getMaxX() - xmax)
        if (bbox.getMaxY() - ymax) < 0:
            ymin -= (bbox.getMaxY() - ymax)
            ymax += (bbox.getMaxY() - ymax)
        if xmin > xmax or ymin > ymax:
            continue
        
        kbbox = afwGeom.Box2I(afwGeom.Point2I(xmin, ymin), afwGeom.Point2I(xmax, ymax))
        try:
            fsb.apply(afwImage.MaskedImageF(templateExposure.getMaskedImage(), kbbox, False).getMask())
            bm1 = fsb.getBits()
            fsb.apply(afwImage.MaskedImageF(scienceExposure.getMaskedImage(), kbbox, False).getMask())
            bm2 = fsb.getBits()
        except Exception, e:
            pass
        else:
            if not((bm1 & badBitMask) or (bm2 & badBitMask)):
                candidateOutList.append(afwDetect.Footprint(kbbox))
    log.info("Selected %d / %d sources for KernelCandidacy" % (len(candidateOutList), len(candidateInList)))
    return candidateOutList
    
#######
# DiaSource filters (here for now)
#######

class DipoleDeblender(object):
    """A functor to deblend a source as a dipole, and return a new
       source with deblended footprints.  This necessarily overrides
       some of the functionality from
       meas_algorithms/python/lsst/meas/algorithms/deblend.py since we
       need a single source that contains the blended peaks, not
       multiple children sources.  This directly calls the core
       deblending code measDeblend.deblend (optionally _fit_psf for
       debugging).
    """
    def __init__(self):
        # Set up defaults to send to deblender
        self.psf_chisq_cut1 = self.psf_chisq_cut2 = self.psf_chisq_cut2b = np.inf # always deblend as Psf
        self.log = pexLog.Log(pexLog.Log.getDefaultLog(), 'lsst.ip.diffim.DipoleDeblender', pexLog.Log.INFO)
        self.sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

    def __call__(self, source, exposure):
        fp     = source.getFootprint()
        peaks  = fp.getPeaks()
        peaksF = [pk.getF() for pk in peaks]
        fbb    = fp.getBBox()
        fmask  = afwImage.MaskU(fbb)
        fmask.setXY0(fbb.getMinX(), fbb.getMinY())
        afwDetect.setMaskFromFootprint(fmask, fp, 1)

        psf        = exposure.getPsf()
        width, height = psf.getKernel().getDimensions()
        psfAttr    = measAlg.PsfAttributes(psf, width//2, height//2)
        psfSigPix  = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT)
        psfFwhmPix = psfSigPix * self.sigma2fwhm 
        subimage   = afwImage.ExposureF(exposure, fbb, True)
        cpsf       = measDeblend.CachingPsf(psf)

        # if fewer than 2 peaks, just return a copy of the source
        if len(peaks) < 2:
            return source.getTable().copyRecord(source)

        # make sure you only deblend 2 peaks; take the brighest and faintest
        speaks = [(p.getPeakValue(), p) for p in peaks]
        speaks.sort() 
        dpeaks = [speaks[0][1], speaks[-1][1]]

        # and only set these peaks in the footprint (peaks is mutable)
        peaks.clear()
        for peak in dpeaks:
            peaks.push_back(peak)

        if True:
            # Call top-level deblend task
            fpres = measDeblend.deblend(fp, exposure.getMaskedImage(), psf, psfFwhmPix,
                                        log = self.log,
                                        psf_chisq_cut1 = self.psf_chisq_cut1,
                                        psf_chisq_cut2 = self.psf_chisq_cut2,
                                        psf_chisq_cut2b = self.psf_chisq_cut2b)
        else:
            # Call lower-level _fit_psf task

            # Prepare results structure
            fpres = measDeblend.PerFootprint()
            fpres.peaks = []
            for pki,pk in enumerate(dpeaks):
                pkres = measDeblend.PerPeak()
                pkres.peak = pk
                pkres.pki = pki
                fpres.peaks.append(pkres)
            
            for pki,(pk,pkres,pkF) in enumerate(zip(dpeaks, fpres.peaks, peaksF)):
                self.log.logdebug('Peak %i' % pki)
                measDeblend._fit_psf(fp, fmask, pk, pkF, pkres, fbb, dpeaks, peaksF, self.log, 
                                     cpsf, psfFwhmPix, 
                                     subimage.getMaskedImage().getImage(), 
                                     subimage.getMaskedImage().getVariance(), 
                                     self.psf_chisq_cut1, self.psf_chisq_cut2, self.psf_chisq_cut2b)


        deblendedSource = source.getTable().copyRecord(source)
        deblendedSource.setParent(source.getId())
        peakList        = deblendedSource.getFootprint().getPeaks()
        peakList.clear()

        for i, peak in enumerate(fpres.peaks):
            if peak.psfflux > 0:
                suffix = "pos"
            else:
                suffix = "neg"
            self.log.info("deblended.centroid.dipole.psf.%s %f %f" % (suffix, peak.center[0], peak.center[1]))
            self.log.info("deblended.chi2dof.dipole.%s %f" % (suffix, peak.chisq / peak.dof))
            self.log.info("deblended.flux.dipole.psf.%s %f" % (suffix, peak.psfflux * np.sum(peak.tmimg.getImage().getArray())))
            peakList.push_back(peak.peak)
        return deblendedSource
        

class SourceFlagChecker(object):
    """A functor to check whether a difference image source has any flags set that should cause it to be labeled bad."""
    def __init__(self, sources, badFlags=['flags.pixel.edge', 'flags.pixel.interpolated.center', 'flags.pixel.saturated.center']):
        self.keys = [sources.getSchema().find(name).key for name in badFlags]
        self.keys.append(sources.table.getCentroidFlagKey())

    def __call__(self, source):
        for k in self.keys:
            if source.get(k):
                return False
        return True

class DipoleChecker(object):
    """A functor to check for dipoles in difference image source tables."""
    def __init__(self):
        pass

    def __call__(self, source):
        return self.getSn(source), self.getCentroid(source), self.getOrientation(source)

    def getSn(self, source):
        posflux = source.get("flux.dipole.psf.pos")
        posfluxErr = source.get("flux.dipole.psf.pos.err")
        negflux = source.get("flux.dipole.psf.neg")
        negfluxErr = source.get("flux.dipole.psf.neg.err")

        # Not a dipole!
        if (posflux < 0) is (negflux < 0):
            return 0

        return np.sqrt((posflux/posfluxErr)**2 + (negflux/negfluxErr)**2)

    def getCentroid(self, source):
        negCen = source.get("flux.dipole.psf.neg.centroid")
        posCen = source.get("flux.dipole.psf.pos.centroid")
        if (False in np.isfinite(negCen)) or (False in np.isfinite(posCen)):
            return None
        
        center = afwGeom.Point2D(0.5*(negCen[0]+posCen[0]),
                                 0.5*(negCen[1]+posCen[1]))
        return center

    def getOrientation(self, source):
        negCen = source.get("flux.dipole.psf.neg.centroid")
        posCen = source.get("flux.dipole.psf.pos.centroid")
        if (False in np.isfinite(negCen)) or (False in np.isfinite(posCen)):
            return None

        dx, dy = posCen[0]-negCen[0], posCen[1]-negCen[1]
        angle  = afwGeom.Angle(np.arctan2(dx, dy), afwGeom.radians)
        return angle

    def displayDipoles(self, exposure, sources, frame=1):
        import lsst.afw.display.ds9 as ds9                
        ds9.mtv(exposure, frame=frame)
        with ds9.Buffering():
            for source in sources:
                cen = source.get("flux.dipole.psf.centroid")
                ctype= "green"
                if (False in np.isfinite(cen)):
                    cen = source.getCentroid()
                    ctype = "red"
                    print "NOT DIPOLE: %.1f,%.1f" % (cen.getX(), cen.getY())
                else:
                    print "DIPOLE: %.1f,%.1f %.1f,%.1f" % (cen.getX(), cen.getY(), source.get("flux.dipole.psf.neg"), source.get("flux.dipole.psf.pos"))

                ds9.dot("o", cen.getX(), cen.getY(), size=2, ctype=ctype, frame=frame)

                negCen = source.get("flux.dipole.psf.neg.centroid")
                posCen = source.get("flux.dipole.psf.pos.centroid")
                if (False in np.isfinite(negCen)) or (False in np.isfinite(posCen)):
                    continue

                ds9.line([(negCen.getX(), negCen.getY()),(posCen.getX(), posCen.getY())], ctype="yellow", frame=frame)
            
        

class NbasisEvaluator(object):
    """A functor to evaluate the Bayesian Information Criterion for the number of basis sets going into the kernel fitting"""
    def __init__(self, psfMatchConfig, psfFwhmPixTc, psfFwhmPixTnc):
        self.psfMatchConfig = psfMatchConfig
        self.psfFwhmPixTc = psfFwhmPixTc
        self.psfFwhmPixTnc = psfFwhmPixTnc
        if not self.psfMatchConfig.kernelBasisSet == "alard-lupton":
            raise RuntimeError, "BIC only implemnted for AL (alard lupton) basis"
        
    def __call__(self, kernelCellSet, log):
        d1, d2, d3 = self.psfMatchConfig.alardDegGauss
        bicArray = {}
        for d1i in range(1, d1+1):
            for d2i in range(1, d2+1):
                for d3i in range(1, d3+1):
                    dList = [d1i, d2i, d3i]
                    bicConfig = type(self.psfMatchConfig)(self.psfMatchConfig, alardDegGauss=dList)
                    kList = makeKernelBasisList(bicConfig, self.psfFwhmPixTc, self.psfFwhmPixTnc)
                    k = len(kList)
                    singlekv = diffimLib.BuildSingleKernelVisitorF(kList, pexConfig.makePolicy(bicConfig))
                    singlekv.setSkipBuilt(False)
                    kernelCellSet.visitCandidates(singlekv, bicConfig.nStarPerCell)        

                    for cell in kernelCellSet.getCellList():
                        for cand in cell.begin(False): # False = include bad candidates
                            cand = diffimLib.cast_KernelCandidateF(cand)
                            if cand.getStatus() != afwMath.SpatialCellCandidate.GOOD:
                                continue
                            diffIm = cand.getDifferenceImage(diffimLib.KernelCandidateF.RECENT)
                            bbox = cand.getKernel(diffimLib.KernelCandidateF.RECENT).shrinkBBox(diffIm.getBBox(afwImage.LOCAL))
                            diffIm = type(diffIm)(diffIm, bbox, True)
                            chi2 = diffIm.getImage().getArray()**2 / diffIm.getVariance().getArray()
                            n = chi2.shape[0] * chi2.shape[1]
                            bic = np.sum(chi2) + k * np.log(n)
                            if not bicArray.has_key(cand.getId()):
                                bicArray[cand.getId()] = {}
                            bicArray[cand.getId()][(d1i, d2i, d3i)] = bic

        bestConfigs = []
        candIds = bicArray.keys()
        for candId in candIds:
            cconfig, cvals = bicArray[candId].keys(), bicArray[candId].values()
            idx = np.argsort(cvals)
            bestConfig = cconfig[idx[0]]
            bestConfigs.append(bestConfig)
        
        counter = Counter(bestConfigs).most_common(3)
        log.info("B.I.C. prefers basis complexity %s %d times; %s %d times; %s %d times" % (counter[0][0], counter[0][1],
                                                                                            counter[1][0], counter[1][1],
                                                                                            counter[2][0], counter[2][1]))
        return counter[0][0], counter[1][0], counter[2][0]
