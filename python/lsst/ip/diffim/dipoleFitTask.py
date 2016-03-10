#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
#import sys
from collections import namedtuple
import numpy as np

from scipy.stats import chi2
import lmfit

#import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
#import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
#import lsst.ip.diffim as ipDiffim
#import lsst.meas.modelfit as measMF
#import lsst.afw.display.ds9 as ds9
import lsst.meas.base as measBase
import lsst.pex.exceptions as pexExcept
from lsst.pex.config import Field
import lsst.pex.logging as pexLog

__all__ = ("DipoleFitConfig", "DipoleFitTask", "DipoleFitPlugin")

## Create a namedtuple to hold all of the relevant output from the lmfit results
resultsOutput = namedtuple('resultsOutput',
                           ['psfFitPosCentroidX', 'psfFitPosCentroidY',
                            'psfFitNegCentroidX', 'psfFitNegCentroidY', 'psfFitPosFlux',
                            'psfFitNegFlux', 'psfFitPosFluxSigma', 'psfFitNegFluxSigma',
                            'psfFitCentroidX', 'psfFitCentroidY', 'psfFitOrientation',
                            'psfFitSignaltoNoise', 'psfFitChi2', 'psfFitRedChi2'])

def dipoleFunc(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=None,
               b=None, x1=None, y1=None, **kwargs):
    """
    dipoleFunc(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg)
    Dipole model generator using difference image's psf.
    Psf is is passed as kwargs['psf']
    Other kwargs include 'rel_weight' - set the relative weighting of pre-sub images
    versus the diffim, and
                         'footprint' - the footprint of the dipole source
    """

    psf = kwargs.get('psf')
    rel_weight = kwargs.get('rel_weight')
    fp = kwargs.get('footprint')
    bbox = fp.getBBox()

    if fluxNeg is None:
        fluxNeg = flux

    gradient = None #gradientImage = None  ## TBD: is using an afwImage faster?
    if b is not None and abs(b)+abs(x1)+abs(y1) > 0.:
        y, x = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
        gradient = b + x1 * x + y1 * y
        # gradientImage = afwImage.ImageF(bbox)
        # gradientImage.getArray()[:,:] = gradient

    p_pos = psf.computeImage(afwGeom.Point2D(xcenPos, ycenPos)).convertF()
    p_pos_sum = np.sum(p_pos.getArray())
    p_pos *= (flux/p_pos_sum)
    p_neg = psf.computeImage(afwGeom.Point2D(xcenNeg, ycenNeg)).convertF()
    p_neg_sum = np.sum(p_neg.getArray())
    p_neg *= (fluxNeg/p_neg_sum)

    pos_box = p_pos.getBBox()
    neg_box = p_neg.getBBox()

    ## Clip the PSF image bounding boxes to fall within the footprint bounding box
    pos_box.clip(bbox)
    neg_box.clip(bbox)

    ## Then actually crop the psf images.
    ## Usually not necessary, but if the dipole is near the edge of the image...
    ## Would be nice if we could compare original pos_box with clipped pos_box and
    ##     see if it actually was clipped.
    p_pos = afwImage.ImageF(p_pos, pos_box, afwImage.PARENT)
    p_neg = afwImage.ImageF(p_neg, neg_box, afwImage.PARENT)

    posIm = afwImage.ImageF(bbox)
    tmpSubim = afwImage.ImageF(posIm, pos_box, afwImage.PARENT)
    tmpSubim += p_pos

    negIm = afwImage.ImageF(bbox)
    tmpSubim = afwImage.ImageF(negIm, neg_box, afwImage.PARENT)
    tmpSubim += p_neg

    if gradient is not None:
        # posIm += gradientImage
        # negIm += gradientImage
        posIm.getArray()[:,:] += gradient
        negIm.getArray()[:,:] += gradient

    diffIm = afwImage.ImageF(bbox)
    diffIm += posIm
    diffIm -= negIm

    zout = diffIm.getArray()
    if rel_weight > 0.:
        zout = np.append([zout], [posIm.getArray(), negIm.getArray()], axis=0)

    return zout

def fitDipole(diffim, source, posImage=None, negImage=None, tol=1e-7, rel_weight=0.1,
              fitBgGradient=True, centroidRangeInSigma=5., separateNegParams=True,
              verbose=False, display=False):
    """
    fitDipole()
    """
    ## diffim is the image difference (exposure)
    ## source is a putative dipole source, with a footprint, from a catalog.
    ## separateNegParams --> separate flux (and TBD: gradient) params for negative img.
    ## Otherwise same as posImage
    ## Returns a lmfit.MinimzerResult object

    fp = source.getFootprint()
    box = fp.getBBox()
    subim = diffim.getMaskedImage().Factory(
        diffim.getMaskedImage(), box, afwImage.PARENT)

    z = subim.getArrays()[0] ## allow passing of just the diffim
    weights = subim.getArrays()[2]  ## get the weights (=1/variance)
    if posImage is not None:
        posSubim = posImage.getMaskedImage().Factory(
            posImage.getMaskedImage(), box, afwImage.PARENT)
        negSubim = negImage.getMaskedImage().Factory(
            negImage.getMaskedImage(), box, afwImage.PARENT)
        z = np.append([z], [posSubim.getArrays()[0],
                        negSubim.getArrays()[0]], axis=0)
        weights = np.append([weights], [posSubim.getArrays()[2] * rel_weight,
                        negSubim.getArrays()[2] * rel_weight], axis=0)

    weights = 1. / weights  ## TBD: is there an inplace operator for this?

    psfSigma = diffim.getPsf().computeShape().getDeterminantRadius()

    ## Create the lmfit model (uses scipy 'leastsq' option by default - Levenberg-Marquardt)
    gmod = lmfit.Model(dipoleFunc, verbose=verbose)

    pks = fp.getPeaks()
    cenPos, cenNeg = pks[0].getF(), pks[1].getF()

    ## For close/faint dipoles the starting locs (min/max) might be way off, let's help them a bit.
    ## First assume dipole is not separated by more than 5*psfSigma.
    centroidRange = psfSigma * centroidRangeInSigma

    ## parameter hints/constraints: https://lmfit.github.io/lmfit-py/model.html#model-param-hints-section
    ## might make sense to not use bounds -- see http://lmfit.github.io/lmfit-py/bounds.html
    gmod.set_param_hint('xcenPos', value=cenPos[0], min=cenPos[0]-centroidRange, max=cenPos[0]+centroidRange)
    gmod.set_param_hint('ycenPos', value=cenPos[1], min=cenPos[1]-centroidRange, max=cenPos[1]+centroidRange)
    gmod.set_param_hint('xcenNeg', value=cenNeg[0], min=cenNeg[0]-centroidRange, max=cenNeg[0]+centroidRange)
    gmod.set_param_hint('ycenNeg', value=cenNeg[1], min=cenNeg[1]-centroidRange, max=cenNeg[1]+centroidRange)

    ## Estimate starting flux. This strongly affects runtime performance so we want to make it close.
    ## Subtract a (constant) background (median) to get estimate - this may not be necessary.
    ## Just using the area within the footprint bounding box.
    if posImage is not None:
        #startingFlux = (z[1,:] - np.median(z[1,:])).sum()   ## use the pos. image
        startingFlux = (np.abs(z[0,:]) - np.median(z[0,:])).sum() / 2.  ## use the dipole
    else:
        startingFlux = (np.abs(z) - np.median(z)).sum() / 2.  ## use the dipole for an estimate.

    ## TBD: set max. flux limit?
    gmod.set_param_hint('flux', value=startingFlux, min=0.1) #, max=startingFlux * 2.)

    if separateNegParams:
        ## TBD: set max negative lobe flux limit?
        gmod.set_param_hint('fluxNeg', value=startingFlux, min=0.1) #, max=startingFlux * 2.)
    else:
        gmod.set_param_hint('fluxNeg', value=startingFlux, min=0.1, expr='flux')

    ## Fixed parameters (dont fit for them if there are no pre-sub images or no gradient fit requested):
    varyBgParams = (rel_weight > 0. and fitBgGradient)
    gmod.set_param_hint('b', value=0., vary=varyBgParams)
    gmod.set_param_hint('x1', value=0., vary=varyBgParams)
    gmod.set_param_hint('y1', value=0., vary=varyBgParams)

    ## Compute footprint bounding box as a numpy extent
    extent = (box.getBeginX(), box.getEndX(), box.getBeginY(), box.getEndY())
    in_x = np.array(extent)   # input x coordinate grid

    ## Note that while we can, we don't need to set initial values for params here, since we
    ## set their param_hint's above.
    ## add "method" param to not use 'leastsq' (==levenberg-marquardt), e.g. "method='nelder'"
    result = gmod.fit(z, weights=weights, x=in_x,
                      verbose=verbose,
                      fit_kws={'ftol':tol, 'xtol':tol, 'gtol':tol}, ## see scipy docs
                      psf=diffim.getPsf(),
                      rel_weight=rel_weight,
                      footprint=fp
                     )

    ## Probably never wanted - also this takes a long time (longer than the fit!)
    ## This is how to get confidence intervals out;
    ##    see here: https://lmfit.github.io/lmfit-py/confidence.html and
    ##    here: http://cars9.uchicago.edu/software/python/lmfit/model.html
    if verbose and separateNegParams:
        print result.fit_report(show_correl=False)
        print result.ci_report()

    ## Display images, model fits and residuals (currently uses matplotlib display functions
    if display:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 2.5))
        plt.subplot(1, 3, 1)
        if posImage is not None:
            display2dArray(z[0,:], 'Data', True, extent=extent)
        else:
            display2dArray(z, 'Data', True, extent=extent)
        plt.subplot(1, 3, 2)
        if posImage is not None:
            display2dArray(result.best_fit[0,:], 'Model', True, extent=extent)
        else:
            display2dArray(result.best_fit, 'Model', True, extent=extent)
        plt.subplot(1, 3, 3)
        if posImage is not None:
            display2dArray(z[0,:] - result.best_fit[0,:], 'Residual', True, extent=extent)
        else:
            display2dArray(z - result.best_fit, 'Residual', True, extent=extent)

    return result

def fitDipole_new(exposure, source, posImage=None, negImage=None, tol=1e-7, rel_weight=0.1,
                  return_fitObj=False, fitBgGradient=True, centroidRangeInSigma=5.,
                  separateNegParams=True, verbose=False, display=False):

    result = fitDipole(exposure, source=source, posImage=posImage, negImage=negImage,
                       rel_weight=rel_weight, tol=tol, fitBgGradient=fitBgGradient,
                       centroidRangeInSigma=centroidRangeInSigma,
                       separateNegParams=separateNegParams, verbose=verbose, display=display)

    results = result.best_values

    centroid = ((results['xcenPos']+results['xcenNeg'])/2., (results['ycenPos']+results['ycenNeg'])/2.)
    dx, dy = results['xcenPos'] - results['xcenNeg'], results['ycenPos'] - results['ycenNeg']
    angle = np.arctan2(dy, dx) / np.pi * 180.   ## convert to degrees (should keep as rad?)

    ## TBD - signalToNoise should be flux / variance_within_footprint, not flux / fluxErr.
    fluxVal, fluxErr = result.params['flux'].value, result.params['flux'].stderr
    fluxValNeg, fluxErrNeg = result.params['fluxNeg'].value, result.params['fluxNeg'].stderr
    signalToNoise = np.sqrt((fluxVal/fluxErr)**2 + (fluxValNeg/fluxErrNeg)**2) ## Derived from DipoleAnalysis

    chi2, redchi2 = result.chisqr, result.redchi

    out = resultsOutput(results['xcenPos'], results['ycenPos'], results['xcenNeg'], results['ycenNeg'],
                        fluxVal, -fluxValNeg, fluxErr, fluxErrNeg, centroid[0], centroid[1], angle,
                        signalToNoise, chi2, redchi2)

    if not return_fitObj:  ## for debugging
        return out
    return out, result

## Create a new measurement task (`DipoleFitTask`) that can handle all other SFM tasks but can
## pass a separate pos- and neg- exposure/image to the `DipoleFitPlugin`s `run()` method.

class DipoleFitConfig(measBase.SingleFramePluginConfig):

    centroidRange = Field(
        dtype=float, default=5.,
        doc="assume dipole is not separated by more than centroidRange*psfSigma")

    relWeight = Field(
        dtype=float, default=1.,
        doc="relative weighting of pre-subtraction images")

    tolerance = Field(
        dtype=float, default=1e-7,
        doc="fit tolerance")

    fitBgGradient = Field(
        dtype=bool, default=True,
        doc="fit separate parameters for linear gradient in pre-sub. images")

    fitSeparateNegParams = Field(
        dtype=bool, default=True,
        doc="fit parameters for negative values (flux/gradient) separately from pos.")

    """!Config params for classification of detected diaSources as dipole or not"""
    minSn = Field(
        dtype=float, default=np.sqrt(2) * 5.0,
        doc="Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole")

    maxFluxRatio = Field(
        dtype = float, default = 0.65,
        doc = "Maximum flux ratio in either lobe to be considered a dipole")

    ## Choose a maxChi2DoF corresponding to a significance level of at most 0.05
    ## (note this is actually a significance not a chi2 number)
    maxChi2DoF = Field(
        dtype = float, default = 0.05,
        doc = "Maximum Chi2/DoF of fit to be considered a dipole")

    verbose = Field(
        dtype=bool, default=False,
        doc="be verbose; this is slow")

class DipoleFitTask(measBase.SingleFrameMeasurementTask):

    ConfigClass = DipoleFitConfig
    _DefaultName = "ip_diffim_DipoleFit"

    def __init__(self, schema, algMetadata=None, **kwds):

        measBase.SingleFrameMeasurementTask.__init__(self, schema, algMetadata, **kwds)

        self.dpFitConfig = DipoleFitConfig()
        self.dipoleFitter = DipoleFitPlugin(self.dpFitConfig, name=self._DefaultName,
                                            schema=schema, metadata=algMetadata)

    def run(self, sources, exposure, posImage=None, negImage=None, **kwds):
        """!Run dipole measurement and classification
        @param sources       diaSources that will be measured using dipole measurement
        @param exposure      Exposure on which the diaSources were detected
        @param **kwds        Sent to SingleFrameMeasurementTask
        """

        measBase.SingleFrameMeasurementTask.run(self, sources, exposure, **kwds)
        #self.dipoleFitter.posImage = posImage
        #self.dipoleFitter.negImage = negImage

        if not sources:
            return

        for source in sources:
            self.dipoleFitter.measure(source, exposure, posImage, negImage)

class DipoleFitTransform(measBase.FluxTransform):

    def __init__(self, config, name, mapper):
        measBase.FluxTransform.__init__(self, name, mapper)
        mapper.addMapping(mapper.getInputSchema().find(name + "_flag_edge").key)

@measBase.register("ip_diffim_DipoleFit")
class DipoleFitPlugin(measBase.SingleFramePlugin):

    ConfigClass = DipoleFitConfig

    FAILURE_EDGE = 1

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER    ## algorithms that require both getShape() and getCentroid(),
                                 ## in addition to a Footprint and its Peaks

    @classmethod
    def getTransformClass(cls):
        return DipoleFitTransform

    def __init__(self, config, name, schema, metadata):
        measBase.SingleFramePlugin.__init__(self, config, name, schema, metadata)

        self.log = pexLog.Log(pexLog.Log.getDefaultLog(),
                              'lsst.ip.diffim.DipoleFitPlugin', pexLog.Log.INFO)

        self._setupSchema(config, name, schema, metadata)

    def _setupSchema(self, config, name, schema, metadata):
        # Get a FunctorKey that can quickly look up the "blessed" centroid value.
        self.centroidKey = afwTable.Point2DKey(schema["slot_Centroid"])

        # Add some fields for our outputs, and save their Keys.
        # Use setAttr() to programmatically set the pos/neg named attributes to values, e.g.
        # self.posCentroidKeyX = 'ip_diffim_DipoleFit_pos_centroid_x'

        for pos_neg in ['pos', 'neg']:

            key = schema.addField(
                schema.join(name, pos_neg, "flux"), type=float, units="dn",
                doc="Dipole {0} lobe flux".format(pos_neg))
            setattr(self, ''.join((pos_neg, 'FluxKey')), key)

            key = schema.addField(
                schema.join(name, pos_neg, "fluxSigma"), type=float, units="pixels",
                doc="1-sigma uncertainty for {0} dipole flux".format(pos_neg))
            setattr(self, ''.join((pos_neg, 'FluxSigmaKey')), key)

            for x_y in ['x', 'y']:
                key = schema.addField(
                    schema.join(name, pos_neg, "centroid", x_y), type=float, units="pixels",
                    doc="Dipole {0} lobe centroid".format(pos_neg))
                setattr(self, ''.join((pos_neg, 'CentroidKey', x_y.upper())), key)

        for x_y in ['x', 'y']:
            key = schema.addField(
                schema.join(name, "centroid", x_y), type=float, units="pixels",
                doc="Dipole centroid")
            setattr(self, ''.join(('centroidKey', x_y.upper())), key)

        self.fluxKey = schema.addField(
                schema.join(name, "flux"), type=float, units="dn",
                doc="Dipole overall flux")

        self.orientationKey = schema.addField(
            schema.join(name, "orientation"), type=float, units="deg",
            doc="Dipole orientation")

        self.chi2dofKey = schema.addField(
            schema.join(name, "chi2dof"), type=float,
            doc="Chi2 per degree of freedom of dipole fit")

        self.classificationFlagKey = schema.addField(
            schema.join(name, "flag", "classification"), type="Flag",
            doc="flag indicating source is classified as being a dipole")

        self.flagKey = schema.addField(
            schema.join(name, "flag"), type="Flag",
            doc="general failure flag for dipole fit")

        self.edgeFlagKey = schema.addField(
            schema.join(name, "flag", "edge"), type="Flag",
            doc="flag set when rectangle used by dipole doesn't fit in the image")

    def measure(self, measRecord, exposure, posImage=None, negImage=None):
        ## Do the non-linear least squares estimation
        try:
            result = fitDipole_new(exposure, measRecord,
                                   posImage=posImage, negImage=negImage,
                                   rel_weight=self.config.relWeight,
                                   tol=self.config.tolerance,
                                   centroidRangeInSigma=self.config.centroidRange,
                                   fitBgGradient=self.config.fitBgGradient,
                                   separateNegParams=self.config.fitSeparateNegParams,
                                   verbose=self.config.verbose, display=False)
        except pexExcept.LengthError as err:
            raise measBase.MeasurementError(err, self.FAILURE_EDGE)

        ## TBD: add chi2, coord/flux uncertainties, dipole classification

        self.log.log(self.log.DEBUG, "Dipole fit result: %s" % str(result))

        ## Add the relevant values to the measRecord
        measRecord[self.posFluxKey] = result.psfFitPosFlux
        measRecord[self.posFluxSigmaKey] = result.psfFitSignaltoNoise   ## to be changed to actual sigma
        measRecord[self.posCentroidKeyX] = result.psfFitPosCentroidX
        measRecord[self.posCentroidKeyY] = result.psfFitPosCentroidY

        measRecord[self.negFluxKey] = result.psfFitNegFlux
        measRecord[self.negFluxSigmaKey] = result.psfFitSignaltoNoise
        measRecord[self.negCentroidKeyX] = result.psfFitNegCentroidX
        measRecord[self.negCentroidKeyY] = result.psfFitNegCentroidY

        ## Dia source flux: average of pos+neg
        measRecord[self.fluxKey] = (abs(result.psfFitPosFlux) + abs(result.psfFitNegFlux))/2.
        measRecord[self.orientationKey] = result.psfFitOrientation
        measRecord[self.centroidKeyX] = (result.psfFitPosCentroidX + result.psfFitNegCentroidX)/2.
        measRecord[self.centroidKeyY] = (result.psfFitPosCentroidY + result.psfFitNegCentroidY)/2.

        measRecord[self.chi2dofKey] = result.psfFitRedChi2

        ## Determine if source is classified as dipole (similar to orig. dipole classification task)
        ## First, does the total signal-to-noise surpass the minSn?
        passesSn = result.psfFitSignaltoNoise > self.config.minSn
        ## Second, does are the pos/neg fluxes no more than 0.65 of the total flux?
        passesFluxPos = abs(result.psfFitPosFlux) / (measRecord[self.fluxKey]*2.) < self.config.maxFluxRatio
        passesFluxNeg = abs(result.psfFitNegFlux) / (measRecord[self.fluxKey]*2.) < self.config.maxFluxRatio

        ## Third, is it a good fit (chi2dof < 1)?
        ## Use scipy's chi2 cumulative distrib to estimate significance
        ndof = result.psfFitChi2 / result.psfFitRedChi2
        significance = chi2.cdf(result.psfFitChi2, ndof)
        passesChi2 = significance < self.config.maxChi2DoF
        #print ndof, result.psfFitChi2, result.psfFitRedChi2, significance, self.config.maxChi2DoF

        if (passesSn and passesFluxPos and passesFluxNeg and passesChi2):
            measRecord.set(self.classificationFlagKey, True)
        else:
            measRecord.set(self.classificationFlagKey, False)

    def fail(self, measRecord, error=None):
        measRecord.set(self.flagKey, True)
        if error is not None:
            assert error.getFlagBit() == self.FAILURE_EDGE
            measRecord.set(self.edgeFlagKey, True)




################# UTILITIES FUNCTIONS -- TBD WHERE THEY ULTIMATELY END UP #######

def display2dArray(arr, title='Data', showBars=True, extent=None):
    import matplotlib.pyplot as plt
    img = plt.imshow(arr, origin='lower', interpolation='none', cmap='gray', extent=extent)
    plt.title(title)
    if showBars:
        plt.colorbar(img, cmap='gray')

def displayImage(image, showBars=True, width=8, height=2.5):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(width, height))
    bbox = image.getBBox()
    extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
    plt.subplot(1, 3, 1)
    ma = image.getArray()
    display2dArray(ma, title='Data', showBars=showBars, extent=extent)

def displayMaskedImage(maskedImage, showMasks=True, showVariance=False, showBars=True, width=8, height=2.5):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(width, height))
    bbox = maskedImage.getBBox()
    extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
    plt.subplot(1, 3, 1)
    ma = maskedImage.getArrays()
    display2dArray(ma[0], title='Data', showBars=showBars, extent=extent)
    if showMasks:
        plt.subplot(1, 3, 2)
        display2dArray(ma[1], title='Masks', showBars=showBars, extent=extent)
    if showVariance:
        plt.subplot(1, 3, 3)
        display2dArray(ma[2], title='Variance', showBars=showBars, extent=extent)

def displayExposure(exposure, showMasks=True, showVariance=False, showPsf=False, showBars=True,
                    width=8, height=2.5):
    displayMaskedImage(exposure.getMaskedImage(), showMasks, showVariance=not showPsf, showBars=showBars,
                       width=width, height=height)
    if showPsf:
        plt.subplot(1, 3, 3)
        psfIm = exposure.getPsf().computeImage()
        bbox = psfIm.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
        display2dArray(psfIm.getArray(), title='PSF', showBars=showBars, extent=extent)

def makeStarImage_lsst(w=101, h=101, xc=[15.3], yc=[18.6], flux=[2500], psfSigma=2., noise=1.0,
                       gradientParams=None, schema=None):
    from lsst.meas.base.tests import TestDataset
    bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Point2I(w-1, h-1))
    dataset = TestDataset(bbox, psfSigma=psfSigma, threshold=1.)
    for i in xrange(len(xc)):
        dataset.addSource(flux=flux[i], centroid=afwGeom.Point2D(xc[i], yc[i]))
    if schema is None:
        schema = TestDataset.makeMinimalSchema()
    exposure, catalog = dataset.realize(noise=noise, schema=schema)
    if gradientParams is not None:
        imgArr = exposure.getMaskedImage().getArrays()[0]
        y, x = np.mgrid[:w, :h]
        gp = gradientParams
        gradient = gp[0] + gp[1] * x + gp[2] * y
        imgArr += gradient

    return exposure, catalog

def makeDipoleImage_lsst(w=101, h=101, xcenPos=[27.], ycenPos=[25.], xcenNeg=[23.], ycenNeg=[25.],
                         psfSigma=2., flux=[30000.], fluxNeg=[None], noise=10., gradientParams=None):

    posImage, posCatalog = makeStarImage_lsst(w, h, xcenPos, ycenPos, flux=flux, psfSigma=psfSigma,
                                              gradientParams=gradientParams, noise=noise)
    if fluxNeg is None:
        fluxNeg = flux
    negImage, negCatalog = makeStarImage_lsst(w, h, xcenNeg, ycenNeg, flux=fluxNeg, psfSigma=psfSigma,
                                              gradientParams=gradientParams, noise=noise)

    dipole = posImage.clone()
    di = dipole.getMaskedImage()
    di -= negImage.getMaskedImage()

    ## Carry through pos/neg detection masks to new planes in diffim image
    dm = di.getMask()
    posDetectedBits = posImage.getMaskedImage().getMask().getArray() == dm.getPlaneBitMask("DETECTED")
    negDetectedBits = negImage.getMaskedImage().getMask().getArray() == dm.getPlaneBitMask("DETECTED")
    pos_det = dm.addMaskPlane("DETECTED_POS") ## new mask plane -- this is different from "DETECTED"
    neg_det = dm.addMaskPlane("DETECTED_NEG") ## new mask plane -- this is different from "DETECTED_NEGATIVE"
    dma = dm.getArray()
    ## set the two custom mask planes to these new masks
    dma[:,:] = posDetectedBits * pos_det + negDetectedBits * neg_det
    return dipole, (posImage, posCatalog), (negImage, negCatalog)

def detectDipoleSources(dipole, posImage, posCatalog, negImage, negCatalog, doMerge=True):
    from lsst.meas.deblender import SourceDeblendTask

    # Start with a minimal schema - only the fields all SourceCatalogs need
    schema = afwTable.SourceTable.makeMinimalSchema()

    # Create and run a task for deblending (optional, but almost always a good idea).
    # Again, the task defines a few flag fields it will later fill.
    deblendTask = SourceDeblendTask(schema=schema)

    deblendTask.run(posImage, posCatalog, psf=posImage.getPsf())
    deblendTask.run(negImage, negCatalog, psf=negImage.getPsf())

    # Customize the detection task a bit (optional)
    detectConfig = measAlg.SourceDetectionConfig()
    detectConfig.returnOriginalFootprints = False # should be the default
    detectConfig.thresholdValue = 10 # only 10-sigma detections

    ## code from imageDifference.py:
    detectConfig.thresholdPolarity = "both"
    detectConfig.thresholdValue = 5.5
    detectConfig.reEstimateBackground = False
    detectConfig.thresholdType = "pixel_stdev"

    # Create the detection task. We pass the schema so the task can declare a few flag fields
    detectTask = measAlg.SourceDetectionTask(config=detectConfig, schema=schema)

    table = afwTable.SourceTable.make(schema)
    detectResult = detectTask.run(table, dipole)
    catalog = detectResult.sources
    #results = detectTask.makeSourceCatalog(table, exposure, sigma=psfSigma)

    deblendTask.run(dipole, catalog, psf=dipole.getPsf())

    ## Now do the merge.
    if doMerge:
        fpSet = detectResult.fpSets.positive
        fpSet.merge(detectResult.fpSets.negative, 2, 2, False)
        sources = afwTable.SourceCatalog(table)
        fpSet.makeSources(sources)

        return sources

    else:
        return detectTask, deblendTask, schema
