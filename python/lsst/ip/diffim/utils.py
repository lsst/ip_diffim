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

"""Support utilities for Measuring sources"""

import re
import sys

import numpy as num

import lsst.pex.exceptions as pexExcept
import lsst.pex.logging as pexLog
import lsst.daf.base as dafBase
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils
import lsst.meas.algorithms as measAlg
from . import diffimLib
from . import diffimTools

keptPlots = False                       # Have we arranged to keep spatial plots open?
    
def showSourceSet(sSet, xy0=(0, 0), frame=0, ctype=ds9.GREEN, symb="+", size=2):
    """Draw the (XAstrom, YAstrom) positions of a set of Sources.  Image has the given XY0"""

    with ds9.Buffering():
        for s in sSet:
            xc, yc = s.getXAstrom() - xy0[0], s.getYAstrom() - xy0[1]

            if symb == "id":
                ds9.dot(str(s.getId()), xc, yc, frame=frame, ctype=ctype, size=size)
            else:
                ds9.dot(symb, xc, yc, frame=frame, ctype=ctype, size=size)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# Kernel display utilities
#
def showKernelSpatialCells(maskedIm, kernelCellSet, showChi2=False, symb="o",
                           ctype=None, ctypeUnused=None, ctypeBad=None, size=3, frame=None):
    """Show the SpatialCells.  If symb is something that ds9.dot understands (e.g. "o"), the top nMaxPerCell candidates will be indicated with that symbol, using ctype and size"""

    ds9.mtv(maskedIm, frame=frame)
    with ds9.Buffering():
        origin = [-maskedIm.getX0(), -maskedIm.getY0()]
        for cell in kernelCellSet.getCellList():
            displayUtils.drawBBox(cell.getBBox(), origin=origin, frame=frame)

            i = 0
            goodies = ctypeBad is None
            for cand in cell.begin(goodies):
                cand = diffimLib.cast_KernelCandidateF(cand)
                xc, yc = cand.getXCenter() + origin[0], cand.getYCenter() + origin[1]
                if cand.getStatus() == afwMath.SpatialCellCandidate.BAD:
                    color = ctypeBad
                elif cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    color = ctype
                elif cand.getStatus() == afwMath.SpatialCellCandidate.UNKNOWN:
                    color = ctypeUnused
                else:
                    continue

                if color:
                    ds9.dot(symb, xc, yc, frame=frame, ctype=color, size=size)

                    if showChi2:
                        rchi2 = cand.getChi2()
                        if rchi2 > 1e100:
                            rchi2 = num.nan
                        ds9.dot("%d %.1f" % (cand.getId(), rchi2),
                                xc - size, yc - size - 4, frame=frame, ctype=color, size=size)


def showDiaSources(sources, exposure, isFlagged, isDipole, frame=None):
    """Display Dia Sources
    """
    #
    # Show us the ccandidates
    #
    # Too many mask planes in diffims
    for plane in ("BAD", "CR", "EDGE", "INTERPOlATED", "INTRP", "SAT", "SATURATED"):
        ds9.setMaskPlaneVisibility(plane, False)

    mos = displayUtils.Mosaic()
    for i in range(len(sources)):
        source = sources[i]
        badFlag = isFlagged[i]
        dipoleFlag = isDipole[i]
        bbox = source.getFootprint().getBBox()
        stamp = exposure.Factory(exposure, bbox, True)
        im = displayUtils.Mosaic(gutter=1, background=0, mode="x")
        im.append(stamp.getMaskedImage())
        lab = "%.1f,%.1f:" % (source.getX(), source.getY())
        if badFlag:
            ctype = ds9.RED
            lab += "BAD"
        if dipoleFlag:
            ctype = ds9.YELLOW
            lab += "DIPOLE"
        if not badFlag and not dipoleFlag:
            ctype = ds9.GREEN
            lab += "OK"
        mos.append(im.makeMosaic(), lab, ctype)
        #print source.getId()
    title = "Dia Sources"
    mosaicImage = mos.makeMosaic(frame=frame, title=title)
    return mosaicImage

def showKernelCandidates(kernelCellSet, kernel, background, frame=None, showBadCandidates=True,
                         resids=False, kernels=False):
    """Display the Kernel candidates.
    If kernel is provided include spatial model and residuals;  
    If chi is True, generate a plot of residuals/sqrt(variance), i.e. chi
    """

    #
    # Show us the ccandidates
    #
    mos = displayUtils.Mosaic()
    #
    candidateCenters = []
    candidateCentersBad = []
    candidateIndex = 0
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False): # include bad candidates
            cand = diffimLib.cast_KernelCandidateF(cand)

            rchi2 = cand.getChi2()
            if rchi2 > 1e100:
                rchi2 = num.nan

            if not showBadCandidates and cand.isBad():
                continue

            im_resid = displayUtils.Mosaic(gutter=1, background=-5, mode="x")

            try:
                im = cand.getScienceMaskedImage()
                im = im.Factory(im, True)
                im.setXY0(cand.getScienceMaskedImage().getXY0())
            except:
                continue
            if (not resids and not kernels):
                im_resid.append(im.Factory(im, True))
            try:
                im = cand.getTemplateMaskedImage()
                im = im.Factory(im, True)
                im.setXY0(cand.getTemplateMaskedImage().getXY0())
            except:
                continue
            if (not resids and not kernels):
                im_resid.append(im.Factory(im, True))

            # Difference image with original basis
            resid = cand.getDifferenceImage(diffimLib.KernelCandidateF.ORIG)
            if resids:
                var = resid.getVariance()
                var = var.Factory(var, True)
                num.sqrt(var.getArray(), var.getArray()) # inplace sqrt
                resid = resid.getImage()
                resid /= var
                bbox = kernel.shrinkBBox(resid.getBBox())
                resid = resid.Factory(resid, bbox, True)
            elif kernels:
                kim = cand.getKernelImage(diffimLib.KernelCandidateF.ORIG).convertF()
                resid = kim.Factory(kim, True)
            im_resid.append(resid)

            # residuals using spatial model
            ski  = afwImage.ImageD(kernel.getDimensions())
            kernel.computeImage(ski, False, int(cand.getXCenter()), int(cand.getYCenter()))
            sk   = afwMath.FixedKernel(ski)
            sbg  = 0.0
            if background:
                sbg = background(int(cand.getXCenter()), int(cand.getYCenter()))
            sresid = cand.getDifferenceImage(sk, sbg)
            resid  = sresid
            if resids:
                resid = sresid.getImage()
                resid /= var
                bbox = kernel.shrinkBBox(resid.getBBox())
                resid = resid.Factory(resid, bbox, True)
            elif kernels:
                kim = ski.convertF()
                resid = kim.Factory(kim, True)
            im_resid.append(resid)

            im = im_resid.makeMosaic()

            lab = "%d chi^2 %.1f" % (cand.getId(), rchi2)
            ctype = ds9.RED if cand.isBad() else ds9.GREEN

            mos.append(im, lab, ctype)

            if False and num.isnan(rchi2):
                ds9.mtv(cand.getScienceMaskedImage.getImage(), title="candidate", frame=1)
                print "rating",  cand.getCandidateRating()

            im = cand.getScienceMaskedImage()
            center = (candidateIndex, cand.getXCenter() - im.getX0(), cand.getYCenter() - im.getY0())
            candidateIndex += 1
            if cand.isBad():
                candidateCentersBad.append(center)
            else:
                candidateCenters.append(center)

    if resids:
        title = "chi Diffim"
    elif kernels:
        title = "Kernels"
    else:
        title = "Candidates & residuals"
    mosaicImage = mos.makeMosaic(frame=frame, title=title)

    return mosaicImage

def showKernelBasis(kernel, frame=None):
    """Display a Kernel's basis images
    """
    mos = displayUtils.Mosaic()

    for k in kernel.getKernelList():
        im = afwImage.ImageD(k.getDimensions())
        k.computeImage(im, False)
        mos.append(im)
    mos.makeMosaic(frame=frame, title="Kernel Basis Images")

    return mos

###############

def plotKernelSpatialModel(kernel, kernelCellSet, showBadCandidates=True, numSample=128, keepPlots=True, maxCoeff = 10):
    """Plot the Kernel spatial model."""

    try:
        import numpy as num
        import matplotlib.pyplot as plt
        import matplotlib.colors
    except ImportError, e:
        print "Unable to import numpy and matplotlib: %s" % e
        return

    x0 = kernelCellSet.getBBox().getBeginX()
    y0 = kernelCellSet.getBBox().getBeginY()

    candPos = list()
    candFits = list()
    badPos = list()
    badFits = list()
    candAmps = list()
    badAmps = list()
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False):
            cand = diffimLib.cast_KernelCandidateF(cand)
            if not showBadCandidates and cand.isBad():
                continue
            candCenter = afwGeom.PointD(cand.getXCenter(), cand.getYCenter())
            try:
                im = cand.getTemplateMaskedImage()
            except Exception, e:
                continue

            targetFits = badFits if cand.isBad() else candFits
            targetPos = badPos if cand.isBad() else candPos
            targetAmps = badAmps if cand.isBad() else candAmps

            # compare original and spatial kernel coefficients
            kp0 = num.array(cand.getKernel(diffimLib.KernelCandidateF.ORIG).getKernelParameters())
            amp = cand.getCandidateRating()

            targetFits = badFits if cand.isBad() else candFits
            targetPos = badPos if cand.isBad() else candPos
            targetAmps = badAmps if cand.isBad() else candAmps

            targetFits.append(kp0)
            targetPos.append(candCenter)
            targetAmps.append(amp)

    numCandidates = len(candFits)
    numBasisFuncs = kernel.getNBasisKernels()

    xGood = num.array([pos.getX() for pos in candPos]) - x0
    yGood = num.array([pos.getY() for pos in candPos]) - y0
    zGood = num.array(candFits)
    ampGood = num.array(candAmps)

    xBad = num.array([pos.getX() for pos in badPos]) - x0
    yBad = num.array([pos.getY() for pos in badPos]) - y0
    zBad = num.array(badFits)
    ampBad = num.array(badAmps)
    numBad = len(badPos)

    xRange = num.linspace(0, kernelCellSet.getBBox().getWidth(), num=numSample)
    yRange = num.linspace(0, kernelCellSet.getBBox().getHeight(), num=numSample)

    if maxCoeff:
        maxCoeff = min(maxCoeff, kernel.getNKernelParameters())
    else:
        maxCoeff = kernel.getNKernelParameters()

    for k in range(maxCoeff):
        func = kernel.getSpatialFunction(k)
        dfGood = zGood[:,k] - num.array([func(pos.getX(), pos.getY()) for pos in candPos])
        yMin = dfGood.min()
        yMax = dfGood.max()
        if numBad > 0:
            dfBad = zBad[:,k] - num.array([func(pos.getX(), pos.getY()) for pos in badPos])
            # Can really screw up the range...
            yMin = min([yMin, dfBad.min()])
            yMax = max([yMax, dfBad.max()])
        yMin -= 0.05 * (yMax - yMin)
        yMax += 0.05 * (yMax - yMin)

        fRange = num.ndarray((len(xRange), len(yRange)))
        for j, yVal in enumerate(yRange):
            for i, xVal in enumerate(xRange):
                fRange[j][i] = func(xVal, yVal)

        fig = plt.figure(k)

        fig.clf()
        try:
            fig.canvas._tkcanvas._root().lift() # == Tk's raise, but raise is a python reserved word
        except:                                 # protect against API changes
            pass

        fig.suptitle('Kernel component %d' % k)

        # LL
        ax = fig.add_axes((0.1, 0.05, 0.35, 0.35))
        vmin = fRange.min() # - 0.05 * num.fabs(fRange.min())
        vmax = fRange.max() # + 0.05 * num.fabs(fRange.max())
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        im = ax.imshow(fRange, aspect='auto', norm=norm,
                       extent=[0, kernelCellSet.getBBox().getWidth()-1, 0, kernelCellSet.getBBox().getHeight()-1])
        ax.set_title('Spatial polynomial')
        plt.colorbar(im, orientation='horizontal', ticks=[vmin, vmax])

        # UL
        ax = fig.add_axes((0.1, 0.55, 0.35, 0.35))
        ax.plot(-2.5*num.log10(candAmps), zGood[:,k], 'b+')
        if numBad > 0:
            ax.plot(-2.5*num.log10(badAmps), zBad[:,k], 'r+')
        ax.set_title("Basis Coefficients")
        ax.set_xlabel("Instr mag")
        ax.set_ylabel("Coeff")

        # LR
        ax = fig.add_axes((0.55, 0.05, 0.35, 0.35))
        ax.set_autoscale_on(False)
        ax.set_xbound(lower=0, upper=kernelCellSet.getBBox().getHeight())
        ax.set_ybound(lower=yMin, upper=yMax)
        ax.plot(yGood, dfGood, 'b+')
        if numBad > 0:
            ax.plot(yBad, dfBad, 'r+')
        ax.axhline(0.0)
        ax.set_title('dCoeff (indiv-spatial) vs. y')

        # UR
        ax = fig.add_axes((0.55, 0.55, 0.35, 0.35))
        ax.set_autoscale_on(False)
        ax.set_xbound(lower=0, upper=kernelCellSet.getBBox().getWidth())
        ax.set_ybound(lower=yMin, upper=yMax)
        ax.plot(xGood, dfGood, 'b+')
        if numBad > 0:
            ax.plot(xBad, dfBad, 'r+')
        ax.axhline(0.0)
        ax.set_title('dCoeff (indiv-spatial) vs. x')

        fig.show()

    global keptPlots
    if keepPlots and not keptPlots:
        # Keep plots open when done
        def show():
            print "%s: Please close plots when done." % __name__
            try:
                plt.show()
            except:
                pass
            print "Plots closed, exiting..."
        import atexit
        atexit.register(show)
        keptPlots = True


def showKernelMosaic(bbox, kernel, nx=7, ny=None, frame=None, title=None, showCenter=True, showEllipticity=True):
    """Show a mosaic of Kernel images. 
    """
    mos = displayUtils.Mosaic()

    x0 = bbox.getBeginX()
    y0 = bbox.getBeginY()
    width = bbox.getWidth()
    height = bbox.getHeight()

    if not ny:
        ny = int(nx*float(height)/width + 0.5)
        if not ny:
            ny = 1

    schema = afwTable.SourceTable.makeMinimalSchema()
    control = measAlg.GaussianCentroidControl()
    centroider = measAlg.MeasureSourcesBuilder().addAlgorithm(control).build(schema)
    sdssShape = measAlg.SdssShapeControl()
    shaper = measAlg.MeasureSourcesBuilder().addAlgorithm(sdssShape).build(schema)
    table = afwTable.SourceTable.make(schema)
    table.defineCentroid(control.name)
    table.defineShape(sdssShape.name)

    centers = []
    shapes = []
    for iy in range(ny):
        for ix in range(nx):
            x = int(ix*(width-1)/(nx-1)) + x0
            y = int(iy*(height-1)/(ny-1)) + y0

            im = afwImage.ImageD(kernel.getDimensions())
            ksum = kernel.computeImage(im, False, x, y)
            lab = "Kernel(%d,%d)=%.2f" % (x, y, ksum) if False else ""
            mos.append(im, lab)
    
            exp = afwImage.makeExposure(afwImage.makeMaskedImage(im))
            w, h = im.getWidth(), im.getHeight()
            cen = afwGeom.PointD(w//2, h//2)
            src = table.makeRecord()
            foot = afwDet.Footprint(exp.getBBox())
            src.setFootprint(foot)

            centroider.apply(src, exp, cen)
            centers.append((src.getX(), src.getY()))

            shaper.apply(src, exp, cen)
            shapes.append((src.getIxx(), src.getIxy(), src.getIyy()))
            
    mos.makeMosaic(frame=frame, title=title if title else "Model Kernel", mode=nx)

    if centers and frame is not None:
        i = 0
        with ds9.Buffering():
            for cen, shape in zip(centers, shapes):
                bbox = mos.getBBox(i); i += 1
                xc, yc = cen[0] + bbox.getMinX(),  cen[1] + bbox.getMinY()
                if showCenter:
                    ds9.dot("+", xc, yc,  ctype=ds9.BLUE, frame=frame)

                if showEllipticity:
                    ixx, ixy, iyy = shape
                    ds9.dot("@:%g,%g,%g" % (ixx, ixy, iyy), xc, yc, frame=frame, ctype=ds9.RED)

    return mos

def plotPixelResiduals(exposure, warpedTemplateExposure, diffExposure, kernelCellSet, kernel, background, testSources, config, 
                       origVariance = False, nptsFull = 1e6, keepPlots = True, titleFs=14):
    candidateResids = []
    spatialResids   = []
    nonfitResids    = []

    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(True): # only look at good ones
            # Be sure
            if not (cand.getStatus() == afwMath.SpatialCellCandidate.GOOD):
                continue

            cand    = diffimLib.cast_KernelCandidateF(cand)
            diffim  = cand.getDifferenceImage(diffimLib.KernelCandidateF.ORIG)
            orig    = cand.getScienceMaskedImage()

            ski     = afwImage.ImageD(kernel.getDimensions())
            kernel.computeImage(ski, False, int(cand.getXCenter()), int(cand.getYCenter()))
            sk      = afwMath.FixedKernel(ski)
            sbg     = background(int(cand.getXCenter()), int(cand.getYCenter()))
            sdiffim = cand.getDifferenceImage(sk, sbg)

            # trim edgs due to convolution
            bbox    = kernel.shrinkBBox(diffim.getBBox())
            tdiffim  = diffim.Factory(diffim, bbox)
            torig    = orig.Factory(orig, bbox)
            tsdiffim = sdiffim.Factory(sdiffim, bbox)

            if origVariance:
                candidateResids.append(num.ravel(tdiffim.getImage().getArray() / num.sqrt(torig.getVariance().getArray())))
                spatialResids.append(num.ravel(tsdiffim.getImage().getArray() / num.sqrt(torig.getVariance().getArray())))
            else:
                candidateResids.append(num.ravel(tdiffim.getImage().getArray() / num.sqrt(tdiffim.getVariance().getArray())))
                spatialResids.append(num.ravel(tsdiffim.getImage().getArray() / num.sqrt(tsdiffim.getVariance().getArray())))

    fullIm   = diffExposure.getMaskedImage().getImage().getArray()
    fullMask = diffExposure.getMaskedImage().getMask().getArray()
    if origVariance:
        fullVar  = exposure.getMaskedImage().getVariance().getArray()
    else:
        fullVar  = diffExposure.getMaskedImage().getVariance().getArray()

    bitmaskBad  = 0
    bitmaskBad |= afwImage.MaskU.getPlaneBitMask('EDGE')
    bitmaskBad |= afwImage.MaskU.getPlaneBitMask('SAT')
    idx = num.where((fullMask & bitmaskBad) == 0)
    stride = int(len(idx[0]) // nptsFull)
    sidx = idx[0][::stride], idx[1][::stride]
    allResids = fullIm[sidx] / num.sqrt(fullVar[sidx])

    testFootprints = diffimTools.sourceToFootprintList(testSources, warpedTemplateExposure, exposure, config, pexLog.getDefaultLog())
    for fp in testFootprints:
        subexp = diffExposure.Factory(diffExposure, fp.getBBox())
        subim  = subexp.getMaskedImage().getImage()
        if origVariance:
            subvar = afwImage.ExposureF(exposure, fp.getBBox()).getMaskedImage().getVariance()
        else:
            subvar = subexp.getMaskedImage().getVariance()
        nonfitResids.append(num.ravel(subim.getArray() / num.sqrt(subvar.getArray())))

    candidateResids = num.ravel(num.array(candidateResids))
    spatialResids   = num.ravel(num.array(spatialResids))
    nonfitResids    = num.ravel(num.array(nonfitResids))

    try:
        import pylab
        from matplotlib.font_manager import FontProperties
    except ImportError, e:
        print "Unable to import pylab: %s" % e
        return
    
    fig = pylab.figure()
    fig.clf()
    try:
        fig.canvas._tkcanvas._root().lift() # == Tk's raise, but raise is a python reserved word
    except:                                 # protect against API changes
        pass
    if origVariance:
        fig.suptitle("Diffim residuals: Normalized by sqrt(input variance)", fontsize=titleFs)
    else:
        fig.suptitle("Diffim residuals: Normalized by sqrt(diffim variance)", fontsize=titleFs)

    sp1 = pylab.subplot(221)
    sp2 = pylab.subplot(222, sharex=sp1, sharey=sp1)
    sp3 = pylab.subplot(223, sharex=sp1, sharey=sp1)
    sp4 = pylab.subplot(224, sharex=sp1, sharey=sp1)
    xs  = num.arange(-5, 5.05, 0.1)
    ys  = 1. / num.sqrt(2 * num.pi) * num.exp( -0.5 * xs**2 )

    sp1.hist(candidateResids, bins=xs, normed=True, alpha=0.5, label="N(%.2f, %.2f)" % (num.mean(candidateResids), num.var(candidateResids)))
    sp1.plot(xs, ys, "r-", lw=2, label="N(0,1)")
    sp1.set_title("Candidates: basis fit", fontsize=titleFs-2)
    sp1.legend(loc=1, fancybox=True, shadow=True, prop = FontProperties(size=titleFs-6))

    sp2.hist(spatialResids, bins=xs, normed=True, alpha=0.5, label="N(%.2f, %.2f)" % (num.mean(spatialResids), num.var(spatialResids)))
    sp2.plot(xs, ys, "r-", lw=2, label="N(0,1)")
    sp2.set_title("Candidates: spatial fit", fontsize=titleFs-2)
    sp2.legend(loc=1, fancybox=True, shadow=True, prop = FontProperties(size=titleFs-6))

    sp3.hist(nonfitResids, bins=xs, normed=True, alpha=0.5, label="N(%.2f, %.2f)" % (num.mean(nonfitResids), num.var(nonfitResids)))
    sp3.plot(xs, ys, "r-", lw=2, label="N(0,1)")
    sp3.set_title("Control sample: spatial fit", fontsize=titleFs-2)
    sp3.legend(loc=1, fancybox=True, shadow=True, prop = FontProperties(size=titleFs-6))

    sp4.hist(allResids, bins=xs, normed=True, alpha=0.5, label="N(%.2f, %.2f)" % (num.mean(allResids), num.var(allResids)))
    sp4.plot(xs, ys, "r-", lw=2, label="N(0,1)")
    sp4.set_title("Full image (subsampled)", fontsize=titleFs-2)
    sp4.legend(loc=1, fancybox=True, shadow=True, prop = FontProperties(size=titleFs-6))

    pylab.setp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize=titleFs-4)
    pylab.setp(sp2.get_xticklabels()+sp2.get_yticklabels(), fontsize=titleFs-4)
    pylab.setp(sp3.get_xticklabels()+sp3.get_yticklabels(), fontsize=titleFs-4)
    pylab.setp(sp4.get_xticklabels()+sp4.get_yticklabels(), fontsize=titleFs-4)

    sp1.set_xlim(-5, 5)
    sp1.set_ylim(0, 0.5)
    fig.show()

    global keptPlots
    if keepPlots and not keptPlots:
        # Keep plots open when done
        def show():
            print "%s: Please close plots when done." % __name__
            try:
                pylab.show()
            except:
                pass
            print "Plots closed, exiting..."
        import atexit
        atexit.register(show)
        keptPlots = True
