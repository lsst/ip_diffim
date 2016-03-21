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

from collections import namedtuple
import numpy as np

## LSST imports
import lsst.meas.base as meas_base

## Only import what's necessary
from lsst.afw.geom import (Box2I, Point2I, Point2D)
from lsst.afw.image import (ImageF, MaskedImageF, PARENT)
from lsst.afw.table import (Point2DKey)
from lsst.pex.exceptions import LengthError
from lsst.pex.logging import Log
from lsst.pex.config import Field
import lsst.afw.detection as afw_det

__all__ = ("DipoleFitConfig", "DipoleFitTask", "DipoleFitPlugin",
           "DipoleFitAlgorithm", "DipolePlotUtils")

## Create a new measurement task (`DipoleFitTask`) that can handle all other SFM tasks but can
## pass a separate pos- and neg- exposure/image to the `DipoleFitPlugin`s `run()` method.


class DipoleFitConfig(meas_base.SingleFramePluginConfig):

    centroidRange = Field(
        dtype=float, default=5.,
        doc="assume dipole is not separated by more than centroidRange*psfSigma")

    relWeight = Field(
        dtype=float, default=0.5,
        doc="relative weighting of pre-subtraction images")

    tolerance = Field(
        dtype=float, default=1e-7,
        doc="fit tolerance")

    fitBgGradient = Field(
        dtype=bool, default=True,
        doc="fit parameters for linear gradient in pre-sub. images")

    fitSeparateNegParams = Field(
        dtype=bool, default=False,
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


class DipoleFitTask(meas_base.SingleFrameMeasurementTask):

    ConfigClass = DipoleFitConfig
    _DefaultName = "ip_diffim_DipoleFit"

    def __init__(self, schema, algMetadata=None, **kwds):

        meas_base.SingleFrameMeasurementTask.__init__(self, schema, algMetadata, **kwds)

        self.dpFitConfig = DipoleFitConfig()
        self.dipoleFitter = DipoleFitPlugin(self.dpFitConfig, name=self._DefaultName,
                                            schema=schema, metadata=algMetadata)

    def run(self, sources, exposure, posImage=None, negImage=None, **kwds):
        """!Run dipole measurement and classification
        @param sources       diaSources that will be measured using dipole measurement
        @param exposure      Exposure on which the diaSources were detected
        @param **kwds        Sent to SingleFrameMeasurementTask
        """

        meas_base.SingleFrameMeasurementTask.run(self, sources, exposure, **kwds)
        #self.dipoleFitter.posImage = posImage
        #self.dipoleFitter.negImage = negImage

        if not sources:
            return

        for source in sources:
            self.dipoleFitter.measure(source, exposure, posImage, negImage)

class DipoleFitTransform(meas_base.FluxTransform):

    def __init__(self, config, name, mapper):
        meas_base.FluxTransform.__init__(self, name, mapper)
        mapper.addMapping(mapper.getInputSchema().find(name + "_flag_edge").key)

@meas_base.register("ip_diffim_DipoleFit")
class DipoleFitPlugin(meas_base.SingleFramePlugin):

    ConfigClass = DipoleFitConfig

    FAILURE_EDGE = 1   ## too close to the edge
    FAILURE_FIT = 2    ## failure in the fitting
    FAILURE_NOT_DIPOLE = 4  ## input source is not a putative dipole to begin with

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER    ## algorithms that require both getShape() and getCentroid(),
                                 ## in addition to a Footprint and its Peaks

    @classmethod
    def getTransformClass(cls):
        return DipoleFitTransform

    def __init__(self, config, name, schema, metadata):
        meas_base.SingleFramePlugin.__init__(self, config, name, schema, metadata)

        self.log = Log(Log.getDefaultLog(), 'lsst.ip.diffim.DipoleFitPlugin', Log.INFO)

        self._setupSchema(config, name, schema, metadata)

    def _setupSchema(self, config, name, schema, metadata):
        # Get a FunctorKey that can quickly look up the "blessed" centroid value.
        self.centroidKey = Point2DKey(schema["slot_Centroid"])

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

        self.separationKey = schema.addField(
            schema.join(name, "separation"), type=float, units="pixels",
            doc="Pixel separation between positive and negative lobes of dipole")

        self.chi2dofKey = schema.addField(
            schema.join(name, "chi2dof"), type=float,
            doc="Chi2 per degree of freedom of dipole fit")

        self.signalToNoiseKey = schema.addField(
            schema.join(name, "signalToNoise"), type=float,
            doc="Estimated signal-to-noise of dipole fit")

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
        pks = measRecord.getFootprint().getPeaks()
        if len(pks) <= 1: ## not a dipole for our analysis
            self.fail(measRecord, meas_base.MeasurementError('not a dipole', self.FAILURE_NOT_DIPOLE))

        ## Do the non-linear least squares estimation
        try:
            result = DipoleFitAlgorithm.fitDipole_new(
                exposure, measRecord,
                posImage=posImage, negImage=negImage,
                rel_weight=self.config.relWeight,
                tol=self.config.tolerance,
                centroidRangeInSigma=self.config.centroidRange,
                fitBgGradient=self.config.fitBgGradient,
                separateNegParams=self.config.fitSeparateNegParams,
                verbose=self.config.verbose, display=False)
        except LengthError as err:
            raise meas_base.MeasurementError(err, self.FAILURE_EDGE)

        self.log.log(self.log.DEBUG, "Dipole fit result: %s" % str(result))

        if result.psfFitPosFlux <= 1.:   ## usually around 0.1 -- the minimun flux allowed -- i.e. bad fit.
            self.fail(measRecord, meas_base.MeasurementError('dipole fit failure', self.FAILURE_FIT))

        ## add chi2, coord/flux uncertainties (TBD), dipole classification

        ## Add the relevant values to the measRecord
        measRecord[self.posFluxKey] = result.psfFitPosFlux
        measRecord[self.posFluxSigmaKey] = result.psfFitSignaltoNoise   ## to be changed to actual sigma!
        measRecord[self.posCentroidKeyX] = result.psfFitPosCentroidX
        measRecord[self.posCentroidKeyY] = result.psfFitPosCentroidY

        measRecord[self.negFluxKey] = result.psfFitNegFlux
        measRecord[self.negFluxSigmaKey] = result.psfFitSignaltoNoise   ## to be changed to actual sigma!
        measRecord[self.negCentroidKeyX] = result.psfFitNegCentroidX
        measRecord[self.negCentroidKeyY] = result.psfFitNegCentroidY

        ## Dia source flux: average of pos+neg
        measRecord[self.fluxKey] = (abs(result.psfFitPosFlux) + abs(result.psfFitNegFlux))/2.
        measRecord[self.orientationKey] = result.psfFitOrientation
        measRecord[self.separationKey] = np.sqrt((result.psfFitPosCentroidX - result.psfFitNegCentroidX)**2. +
                                                 (result.psfFitPosCentroidY - result.psfFitNegCentroidY)**2.)
        measRecord[self.centroidKeyX] = (result.psfFitPosCentroidX + result.psfFitNegCentroidX)/2.
        measRecord[self.centroidKeyY] = (result.psfFitPosCentroidY + result.psfFitNegCentroidY)/2.

        measRecord[self.signalToNoiseKey] = result.psfFitSignaltoNoise
        measRecord[self.chi2dofKey] = result.psfFitRedChi2

        self.doClassify(measRecord, result)

    def doClassify(self, measRecord, result):
        ## Determine if source is classified as dipole (similar to orig. dipole classification task)
        ## First, does the total signal-to-noise surpass the minSn?
        passesSn = measRecord[self.signalToNoiseKey] > self.config.minSn

        ## Second, does are the pos/neg fluxes no more than 0.65 of the total flux?
        ## By default this will never happen since posFlux = negFlux.
        passesFluxPos = (abs(measRecord[self.posFluxKey]) /
                         (measRecord[self.fluxKey]*2.)) < self.config.maxFluxRatio
        passesFluxNeg = (abs(measRecord[self.negFluxKey]) /
                         (measRecord[self.fluxKey]*2.)) < self.config.maxFluxRatio

        ## Third, is it a good fit (chi2dof < 1)?
        ## Use scipy's chi2 cumulative distrib to estimate significance
        ## This doesn't really work since I don't trust the values in the variance plane (which
        ##   affects the least-sq weights, which affects the resulting chi2).
        from scipy.stats import chi2

        ndof = result.psfFitChi2 / measRecord[self.chi2dofKey]
        significance = chi2.cdf(result.psfFitChi2, ndof)
        passesChi2 = significance < self.config.maxChi2DoF

        allPass = (passesSn and passesFluxPos and passesFluxNeg) # and passesChi2)
        if allPass:  ## Note cannot pass `allPass` into the `measRecord.set()` call below...?
            measRecord.set(self.classificationFlagKey, True)
        else:
            measRecord.set(self.classificationFlagKey, False)

    ## TBD: need to catch more exceptions
    def fail(self, measRecord, error=None):
        measRecord.set(self.flagKey, True)
        if error is not None and error.getFlagBit() == self.FAILURE_EDGE:
            measRecord.set(self.edgeFlagKey, True)


class DipoleFitAlgorithm():
    import lmfit  ## In the future, we might need to change fitters. Astropy, or just scipy, or iminuit?

    ## This is just a private version number to sync with the ipython notebooks that I have been
    ## using for algorithm development.
    __private_version__ = '0.0.1b'

    ## Create a namedtuple to hold all of the relevant output from the lmfit results
    resultsOutput = namedtuple('resultsOutput',
                               ['psfFitPosCentroidX', 'psfFitPosCentroidY',
                                'psfFitNegCentroidX', 'psfFitNegCentroidY', 'psfFitPosFlux',
                                'psfFitNegFlux', 'psfFitPosFluxSigma', 'psfFitNegFluxSigma',
                                'psfFitCentroidX', 'psfFitCentroidY', 'psfFitOrientation',
                                'psfFitSignaltoNoise', 'psfFitChi2', 'psfFitRedChi2'])

    @staticmethod
    def genBgGradientModel(bbox, b=None, x1=0., y1=0., xy=None, x2=0., y2=0.):
        gradient = None #gradientImage = None  ## TBD: is using an afwImage faster?
        if b is not None: ## Don't fit for other gradient parameters if the intercept is not allowed.
            y, x = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
            gradient = np.full_like(x, b, dtype='float64')
            if x1 is not None: gradient += x1 * x
            if y1 is not None: gradient += y1 * y
            if xy is not None: gradient += xy * (x * y)
            if x2 is not None: gradient += x2 * (x * x)
            if y2 is not None: gradient += y2 * (y * y)
            # gradientImage = ImageF(bbox)
            # gradientImage.getArray()[:,:] = gradient
        return gradient

    @staticmethod
    def genStarModel(psf, xcen, ycen, flux, fp):
        ## Generate the psf image, normalize to flux
        psf_img = psf.computeImage(Point2D(xcen, ycen)).convertF()
        psf_img_sum = np.nansum(psf_img.getArray())
        psf_img *= (flux/psf_img_sum)

        ## Clip the PSF image bounding box to fall within the footprint bounding box
        bbox = fp.getBBox()
        psf_box = psf_img.getBBox()
        psf_box.clip(bbox)
        psf_img = ImageF(psf_img, psf_box, PARENT)

        ## Then actually crop the psf image.
        ## Usually not necessary, but if the dipole is near the edge of the image...
        ## Would be nice if we could compare original pos_box with clipped pos_box and
        ##     see if it actually was clipped.
        p_Im = ImageF(bbox)
        tmpSubim = ImageF(p_Im, psf_box, PARENT)
        tmpSubim += psf_img

        return p_Im

    @staticmethod
    def genDipoleModel(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=None,
                       b=None, x1=None, y1=None, xy=None, x2=None, y2=None,
                       bNeg=None, x1Neg=None, y1Neg=None, xyNeg=None, x2Neg=None, y2Neg=None,
                       debug=False, **kwargs):
        """
        genDipoleModel(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg)
        Dipole model generator functor using difference image's psf.
        Psf is is passed as kwargs['psf']
        Other kwargs include 'rel_weight' - set the relative weighting of pre-sub images
        versus the diffim, and
                             'footprint' - the footprint of the dipole source
        Output - generate an output model containing a diffim and the two pre-sub. images
        """

        psf = kwargs.get('psf')
        rel_weight = kwargs.get('rel_weight') ## only says we're including pre-sub. images (if > 0)
        fp = kwargs.get('footprint')
        bbox = fp.getBBox()

        if fluxNeg is None:
            fluxNeg = flux

        if debug:
            print '%.2f %.2f %.2f %.2f %.2f %.2f' % (flux, fluxNeg, xcenPos, ycenPos, xcenNeg, ycenNeg)
            if x1 is not None:
                print '     %.2f %.2f %.2f' % (b, x1, y1)
            if xy is not None:
                print '     %.2f %.2f %.2f' % (xy, x2, y2)

        posIm = DipoleFitAlgorithm.genStarModel(psf, xcenPos, ycenPos, flux, fp)
        negIm = DipoleFitAlgorithm.genStarModel(psf, xcenNeg, ycenNeg, fluxNeg, fp)

        gradient = DipoleFitAlgorithm.genBgGradientModel(bbox, b, x1, y1, xy, x2, y2)
        gradientNeg = gradient
        if bNeg is not None: # and abs(bNeg)+abs(x1Neg)+abs(y1Neg) > 0.:
            gradientNeg = DipoleFitAlgorithm.genBgGradientModel(bbox, bNeg, x1Neg, y1Neg, xyNeg, x2Neg, y2Neg)

        if gradient is not None:
            posIm.getArray()[:,:] += gradient
            negIm.getArray()[:,:] += gradientNeg

        ## Generate the diffIm model
        diffIm = ImageF(bbox)
        diffIm += posIm
        diffIm -= negIm

        zout = diffIm.getArray()
        if rel_weight > 0.:
            zout = np.append([zout], [posIm.getArray(), negIm.getArray()], axis=0)

        return zout

    @staticmethod
    def fitDipole(diffim, source, posImage=None, negImage=None, tol=1e-7, rel_weight=0.5,
                  fitBgGradient=True, bgGradientOrder=1, centroidRangeInSigma=5.,
                  separateNegParams=True, verbose=False, display=False):
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
        subim = MaskedImageF(diffim.getMaskedImage(), box, PARENT)

        z = subim.getArrays()[0] ## allow passing of just the diffim
        weights = subim.getArrays()[2]  ## get the weights (=1/variance)
        if posImage is not None and rel_weight > 0.:
            posSubim = MaskedImageF(posImage.getMaskedImage(), box, PARENT)
            negSubim = MaskedImageF(negImage.getMaskedImage(), box, PARENT)
            z = np.append([z], [posSubim.getArrays()[0],
                                negSubim.getArrays()[0]], axis=0)
            weights = np.append([weights], [posSubim.getArrays()[2] * rel_weight,
                                            negSubim.getArrays()[2] * rel_weight], axis=0)

        weights[:] = 1. / weights  ## TBD: is there an inplace operator for this?

        psfSigma = diffim.getPsf().computeShape().getDeterminantRadius()

        ## Create the lmfit model (lmfit uses scipy 'leastsq' option by default - Levenberg-Marquardt)
        gmod = DipoleFitAlgorithm.lmfit.Model(DipoleFitAlgorithm.genDipoleModel, verbose=verbose)

        ## Add the constraints for centroids, fluxes.
        ## starting constraint - near centroid of footprint
        fpCentroid = np.array([fp.getCentroid().getX(), fp.getCentroid().getY()])
        cenNeg = cenPos = fpCentroid

        pks = fp.getPeaks()
        if len(pks) >= 1:
            cenPos = pks[0].getF()    ## if individual (merged) peaks were detected, use those
        if len(pks) >= 2:
            cenNeg = pks[1].getF()

        ## For close/faint dipoles the starting locs (min/max) might be way off, let's help them a bit.
        ## First assume dipole is not separated by more than 5*psfSigma.
        centroidRange = psfSigma * centroidRangeInSigma

        ## Note - this may be a cheat to assume the dipole is centered in center of the footprint.
        if np.sum(np.sqrt((np.array(cenPos) - fpCentroid)**2.)) > centroidRange:
            cenPos = fpCentroid
        if np.sum(np.sqrt((np.array(cenNeg) - fpCentroid)**2.)) > centroidRange:
            cenPos = fpCentroid

        ## parameter hints/constraints: https://lmfit.github.io/lmfit-py/model.html#model-param-hints-section
        ## might make sense to not use bounds -- see http://lmfit.github.io/lmfit-py/bounds.html
        ## also see this discussion -- https://github.com/scipy/scipy/issues/3129
        gmod.set_param_hint('xcenPos', value=cenPos[0],
                            min=cenPos[0]-centroidRange, max=cenPos[0]+centroidRange)
        gmod.set_param_hint('ycenPos', value=cenPos[1],
                            min=cenPos[1]-centroidRange, max=cenPos[1]+centroidRange)
        gmod.set_param_hint('xcenNeg', value=cenNeg[0],
                            min=cenNeg[0]-centroidRange, max=cenNeg[0]+centroidRange)
        gmod.set_param_hint('ycenNeg', value=cenNeg[1],
                            min=cenNeg[1]-centroidRange, max=cenNeg[1]+centroidRange)

        ## Estimate starting flux. This strongly affects runtime performance so we want to make it close.
        ## I tried many possibilities, the best one seems to be just using sum(abs(diffIm))/2.
        ## I am leaving the other ties here (commented out) for posterity and possible future improvements...

        ## Value to convert peak value to total flux based on flux within psf
        # psfImg = diffim.getPsf().computeImage()
        # pkToFlux = np.nansum(psfImg.getArray()) / diffim.getPsf().computePeak()

        bg = np.nanmedian(z[0,:]) ## Compute the dipole background (probably very close to zero)
        # startingPk = np.nanmax(z[0,:]) - bg   ## use the dipole peak for an estimate.
        # posFlux, negFlux = startingPk * pkToFlux, -startingPk * pkToFlux

        # if len(pks) >= 1:
        #     posFlux = pks[0].getPeakValue() * pkToFlux
        # if len(pks) >= 2:
        #     negFlux = pks[1].getPeakValue() * pkToFlux

        startingFlux = np.nansum(np.abs(z[0,:] - bg)) / 2.   ## use the flux under the dipole for an estimate.
        posFlux = negFlux = startingFlux

        # ## This will only be accurate if there is not a bright gradient/background in the pre-sub images
        # if posImage is not None:
        #     posSubim = ImageF(posImage.getMaskedImage().getImage(), box, PARENT)
        #     negSubim = ImageF(negImage.getMaskedImage().getImage(), box, PARENT)
        #     w, h = posSubim.getWidth(), posSubim.getHeight()

        #     ## If the brightest pixel value is close to a corner's pixel value, then it is picking up the bg gradient
        #     ## In that case, use the footprint peak value instead
        #     if len(pks) >= 1:
        #         posArr = posSubim.getArray()
        #         posPk = np.nanmax(posArr)
        #         if (np.max(np.abs(posPk - np.array([posArr[0,0], posArr[h-1,0],
        #                                             posArr[h-1,w-1], posArr[0,w-1]]))) >
        #             np.abs(posPk - pks[0].getPeakValue())):
        #             posFlux = posPk * pkToFlux

        #     if len(pks) >= 2:
        #         negArr = negSubim.getArray()
        #         negPk = np.nanmax(negArr)
        #         if (np.max(np.abs(negPk - np.array([negArr[0,0], negArr[h-1,0],
        #                                             negArr[h-1,w-1], negArr[0,w-1]]))) >
        #             np.abs(negPk - -pks[1].getPeakValue())):
        #             negFlux = -negPk * pkToFlux

        ## TBD: set max. flux limit?
        gmod.set_param_hint('flux', value=posFlux, min=0.1) #, max=posFlux * 2.)

        if separateNegParams:
            ## TBD: set max negative lobe flux limit?
            gmod.set_param_hint('fluxNeg', value=np.abs(negFlux), min=0.1) #, max=negFlux * 2.)

        ## Fixed parameters (dont fit for them if there are no pre-sub images or no gradient fit requested):
        if (rel_weight > 0. and fitBgGradient):
            if bgGradientOrder >= 0:
                gmod.set_param_hint('b', value=0.)
                if separateNegParams:
                    gmod.set_param_hint('bNeg', value=0.)
            if bgGradientOrder >= 1:
                gmod.set_param_hint('x1', value=0.)
                gmod.set_param_hint('y1', value=0.)
                if separateNegParams:
                    gmod.set_param_hint('x1Neg', value=0.)
                    gmod.set_param_hint('y1Neg', value=0.)
            if bgGradientOrder >= 2:
                gmod.set_param_hint('xy', value=0.)
                gmod.set_param_hint('x2', value=0.)
                gmod.set_param_hint('y2', value=0.)
                if separateNegParams:
                    gmod.set_param_hint('xyNeg', value=0.)
                    gmod.set_param_hint('x2Neg', value=0.)
                    gmod.set_param_hint('y2Neg', value=0.)

        ## Compute footprint bounding box as a numpy extent
        extent = (box.getBeginX(), box.getEndX(), box.getBeginY(), box.getEndY())
        in_x = np.array(extent)   # input x coordinate grid

        if posImage is not None and rel_weight > 0.:
            weights = np.array([np.ones_like(z[0,:]), np.ones_like(z[0,:])*rel_weight,
                                np.ones_like(z[0,:])*rel_weight])
        else:
            weights = 1.

        ## Note that although we can, we're not required to set initial values for params here,
        ## since we set their param_hint's above.
        ## add "method" param to not use 'leastsq' (==levenberg-marquardt), e.g. "method='nelder'"
        result = gmod.fit(z, weights=weights, x=in_x,
                          verbose=verbose,
                          fit_kws={'ftol':tol, 'xtol':tol, 'gtol':tol, 'maxfev':250}, ## see scipy docs
                          psf=diffim.getPsf(),   ## hereon: additional kwargs get passed to genDipoleModel()
                          rel_weight=rel_weight,
                          footprint=fp)

        ## Probably never wanted - also this takes a long time (longer than the fit!)
        ## This is how to get confidence intervals out:
        ##    https://lmfit.github.io/lmfit-py/confidence.html and
        ##    http://cars9.uchicago.edu/software/python/lmfit/model.html
        if verbose:  ## fails if neg params are constrained for some reason.
            print result.fit_report(show_correl=False)
            if separateNegParams:
                print result.ci_report()

        ## Display images, model fits and residuals (currently uses matplotlib display functions)
        if display:
            try:
                DipolePlotUtils.displayFitResults(result, fp)
            except Exception as err:
                print 'Uh oh!', err
                pass

        return result

    @staticmethod
    def fitDipole_new(exposure, source, posImage=None, negImage=None, tol=1e-7, rel_weight=0.1,
                      fitBgGradient=True, centroidRangeInSigma=5., separateNegParams=True,
                      bgGradientOrder=1, verbose=False, display=False, return_fitObj=False):

        fitResult = DipoleFitAlgorithm.fitDipole(
            exposure, source=source, posImage=posImage, negImage=negImage,
            tol=tol, rel_weight=rel_weight, fitBgGradient=fitBgGradient,
            centroidRangeInSigma=centroidRangeInSigma, separateNegParams=separateNegParams,
            bgGradientOrder=bgGradientOrder, verbose=verbose, display=display)

        ## In (rare) extreme cases of very faint dipoles on top of a very steep background
        ## gradient (and we're including the pre-sub. images in the fit), the fit can go haywire.
        ## In this case, hope that the diffim is better and let's re-run the fit just using just
        ## the diffim instead.
        ## This will be rare and running it without background gradient fitting on is about 2x faster
        ## so doesn't add too much time to the fitting.
        ## Currently not run:
        if False and rel_weight > 0. and (fitResult.redchi > 10000. or
                                          fitResult.params['flux'].stderr == 0. or
                                          fitResult.params['flux'].stderr >= 1e6 or
                                          (separateNegParams and fitResult.params['fluxNeg'].stderr == 0)):
            fitResult = DipoleFitAlgorithm.fitDipole(
                exposure, source=source, posImage=posImage, negImage=negImage,
                tol=tol, rel_weight=0., fitBgGradient=False,
                centroidRangeInSigma=centroidRangeInSigma, separateNegParams=separateNegParams,
                bgGradientOrder=0, verbose=verbose, display=display)

        fitParams = fitResult.best_values
        if fitParams['flux'] <= 1.:   ## usually around 0.1 -- the minimun flux allowed -- i.e. bad fit.
            out = DipoleFitAlgorithm.resultsOutput(
                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                np.nan, np.nan, np.nan)
            if return_fitObj:  ## for debugging
                return out, fitResult
            return out

        centroid = ((fitParams['xcenPos']+fitParams['xcenNeg'])/2.,
                    (fitParams['ycenPos']+fitParams['ycenNeg'])/2.)
        dx, dy = fitParams['xcenPos'] - fitParams['xcenNeg'], fitParams['ycenPos'] - fitParams['ycenNeg']
        angle = np.arctan2(dy, dx) / np.pi * 180.   ## convert to degrees (should keep as rad?)

        ## TBD - signalToNoise should be flux / variance_within_footprint, not flux / fluxErr.
        fluxVal, fluxErr = fitParams['flux'], fitResult.params['flux'].stderr
        try:
            fluxValNeg, fluxErrNeg = fitParams['fluxNeg'], fitResult.params['fluxNeg'].stderr
        except:
            fluxValNeg, fluxErrNeg = fitParams['flux'], fitResult.params['flux'].stderr
        try:
            signalToNoise = np.sqrt((fluxVal/fluxErr)**2 + (fluxValNeg/fluxErrNeg)**2) ## Derived from DipoleAnalysis
        except:
            signalToNoise = np.nan

        out = DipoleFitAlgorithm.resultsOutput(
            fitParams['xcenPos'], fitParams['ycenPos'], fitParams['xcenNeg'], fitParams['ycenNeg'],
            fluxVal, -fluxValNeg, fluxErr, fluxErrNeg, centroid[0], centroid[1], angle,
            signalToNoise, fitResult.chisqr, fitResult.redchi)

        if return_fitObj:  ## for debugging
            return out, fitResult
        return out


################# UTILITIES FUNCTIONS -- TBD WHERE THEY ULTIMATELY END UP #######

class DipolePlotUtils():
    try:
        import matplotlib.pyplot as plt
    except Exception as err:
        print 'Uh oh! need matplotlib to use these funcs', err
        pass  ## matplotlib not installed -- cannot do any plotting

    @staticmethod
    def display2dArray(arr, title='Data', showBars=True, extent=None):
        fig = DipolePlotUtils.plt.imshow(arr, origin='lower', interpolation='none', cmap='gray', extent=extent)
        DipolePlotUtils.plt.title(title)
        if showBars:
            DipolePlotUtils.plt.colorbar(fig, cmap='gray')
        return fig

    @staticmethod
    def displayImage(image, showBars=True, width=8, height=2.5):
        fig = DipolePlotUtils.plt.figure(figsize=(width, height))
        bbox = image.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
        DipolePlotUtils.plt.subplot(1, 3, 1)
        ma = image.getArray()
        DipolePlotUtils.display2dArray(ma, title='Data', showBars=showBars, extent=extent)
        return fig

    @staticmethod
    def displayImages(images, showBars=True, width=8, height=2.5):
        fig = DipolePlotUtils.plt.figure(figsize=(width, height))
        for i,image in enumerate(images):
            bbox = image.getBBox()
            extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
            DipolePlotUtils.plt.subplot(1, len(images), i+1)
            ma = image.getArray()
            DipolePlotUtils.display2dArray(ma, title='Data', showBars=showBars, extent=extent)
        return fig

    @staticmethod
    def displayMaskedImage(maskedImage, showMasks=True, showVariance=False, showBars=True, width=8, height=2.5):
        fig = DipolePlotUtils.plt.figure(figsize=(width, height))
        bbox = maskedImage.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
        DipolePlotUtils.plt.subplot(1, 3, 1)
        ma = maskedImage.getArrays()
        DipolePlotUtils.display2dArray(ma[0], title='Data', showBars=showBars, extent=extent)
        if showMasks:
            DipolePlotUtils.plt.subplot(1, 3, 2)
            DipolePlotUtils.display2dArray(ma[1], title='Masks', showBars=showBars, extent=extent)
        if showVariance:
            DipolePlotUtils.plt.subplot(1, 3, 3)
            DipolePlotUtils.display2dArray(ma[2], title='Variance', showBars=showBars, extent=extent)
        return fig

    @staticmethod
    def displayExposure(exposure, showMasks=True, showVariance=False, showPsf=False, showBars=True,
                        width=8, height=2.5):
        fig = DipolePlotUtils.displayMaskedImage(exposure.getMaskedImage(), showMasks,
                                                 showVariance=not showPsf,
                                                 showBars=showBars, width=width, height=height)
        if showPsf:
            DipolePlotUtils.plt.subplot(1, 3, 3)
            psfIm = exposure.getPsf().computeImage()
            bbox = psfIm.getBBox()
            extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
            DipolePlotUtils.display2dArray(psfIm.getArray(), title='PSF', showBars=showBars, extent=extent)
        return fig

    @staticmethod
    def displayCutouts(source, exposure, posImage=None, negImage=None, asHeavyFootprint=False):
        fp = source.getFootprint()
        bbox = fp.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())

        fig = DipolePlotUtils.plt.figure(figsize=(8, 2.5))
        if not asHeavyFootprint:
            subexp = ImageF(exposure.getMaskedImage().getImage(), bbox, PARENT)
        else:
            hfp = afw_det.HeavyFootprintF(fp, exposure.getMaskedImage())
            subexp = DipolePlotUtils.getHeavyFootprintSubimage(hfp)
        DipolePlotUtils.plt.subplot(1, 3, 1)
        DipolePlotUtils.display2dArray(subexp.getArray(), title='Diffim', extent=extent)
        if posImage is not None:
            if not asHeavyFootprint:
                subexp = ImageF(posImage.getMaskedImage().getImage(), bbox, PARENT)
            else:
                hfp = afw_det.HeavyFootprintF(fp, posImage.getMaskedImage())
                subexp = DipolePlotUtils.getHeavyFootprintSubimage(hfp)
            DipolePlotUtils.plt.subplot(1, 3, 2)
            DipolePlotUtils.display2dArray(subexp.getArray(), title='Pos', extent=extent)
        if negImage is not None:
            if not asHeavyFootprint:
                subexp = ImageF(negImage.getMaskedImage().getImage(), bbox, PARENT)
            else:
                hfp = afw_det.HeavyFootprintF(fp, negImage.getMaskedImage())
                subexp = DipolePlotUtils.getHeavyFootprintSubimage(hfp)
            DipolePlotUtils.plt.subplot(1, 3, 3)
            DipolePlotUtils.display2dArray(subexp.getArray(), title='Neg', extent=extent)
        return fig

    @staticmethod
    def makeHeavyCatalog(catalog, exposure, verbose=False):
        """Usage: catalog = makeHeavyCatalog(catalog, exposure, True)"""
        for source in catalog:
            fp = source.getFootprint()
            if not fp.isHeavy():
                if verbose:
                    print 'not heavy => heavy'
                hfp = afw_det.HeavyFootprintF(fp, exposure.getMaskedImage())
                source.setFootprint(hfp)

        return catalog

    @staticmethod
    def getHeavyFootprintSubimage(fp, badfill=np.nan):
        hfp = afw_det.HeavyFootprintF_cast(fp)
        bbox = hfp.getBBox()

        subim2 = ImageF(bbox, badfill)  ## set the mask to NA (can use 0. if desired)
        #subim2.getArray()[:,:] = np.nan
        afw_det.expandArray(hfp, hfp.getImageArray(), subim2.getArray(), bbox.getCorners()[0])
        return subim2

    @staticmethod
    def searchCatalog(catalog, x, y):
        for i,s in enumerate(catalog):
            bbox = s.getFootprint().getBBox()
            if bbox.contains(Point2I(x, y)):
                print i
                return s
        return None

    @staticmethod
    def displayFitResults(result, footprint):
        ## Display images, model fits and residuals (currently uses matplotlib display functions)
        z = result.data
        fit = result.best_fit
        bbox = footprint.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
        if z.shape[0] == 3:
            fig = DipolePlotUtils.plt.figure(figsize=(8, 8))
            for i in range(3):
                DipolePlotUtils.plt.subplot(3, 3, i*3+1)
                DipolePlotUtils.display2dArray(z[i,:], 'Data', True, extent=extent)
                DipolePlotUtils.plt.subplot(3, 3, i*3+2)
                DipolePlotUtils.display2dArray(fit[i,:], 'Model', True, extent=extent)
                DipolePlotUtils.plt.subplot(3, 3, i*3+3)
                DipolePlotUtils.display2dArray(z[i,:] - fit[i,:], 'Residual', True, extent=extent)
            return fig
        else:
            fig = DipolePlotUtils.plt.figure(figsize=(8, 2.5))
            DipolePlotUtils.plt.subplot(1, 3, 1)
            DipolePlotUtils.display2dArray(z, 'Data', True, extent=extent)
            DipolePlotUtils.plt.subplot(1, 3, 2)
            DipolePlotUtils.display2dArray(fit, 'Model', True, extent=extent)
            DipolePlotUtils.plt.subplot(1, 3, 3)
            DipolePlotUtils.display2dArray(z - fit, 'Residual', True, extent=extent)
            return fig

