from __future__ import absolute_import, division, print_function
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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

# LSST imports
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.meas.base as measBase
import lsst.afw.table as afwTable
import lsst.pex.exceptions
import lsst.pex.logging
import lsst.pex.config as pexConfig

__all__ = ("DipoleFitConfig", "DipoleFitTask", "DipoleFitPlugin",
           "DipoleFitAlgorithm", "DipolePlotUtils")


# Create a new measurement task (`DipoleFitTask`) that can handle all other SFM tasks but can
# pass a separate pos- and neg- exposure/image to the `DipoleFitPlugin`s `run()` method.


class DipoleFitConfig(measBase.SingleFramePluginConfig):
    """
    Class to initialize and store dipole fitting configuration parameters
    """

    centroidRange = pexConfig.Field(
        dtype=float, default=5.,
        doc="assume dipole is not separated by more than centroidRange*psfSigma")

    relWeight = pexConfig.Field(
        dtype=float, default=0.5,
        doc="relative weighting of pre-subtraction images")

    tolerance = pexConfig.Field(
        dtype=float, default=1e-7,
        doc="fit tolerance")

    fitBgGradient = pexConfig.Field(
        dtype=bool, default=True,
        doc="fit parameters for linear gradient in pre-sub. images")

    fitSeparateNegParams = pexConfig.Field(
        dtype=bool, default=False,
        doc="fit parameters for negative values (flux/gradient) separately from pos.")

    """!Config params for classification of detected diaSources as dipole or not"""
    minSn = pexConfig.Field(
        dtype=float, default=np.sqrt(2) * 5.0,
        doc="Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole")

    maxFluxRatio = pexConfig.Field(
        dtype=float, default=0.65,
        doc="Maximum flux ratio in either lobe to be considered a dipole")

    # Choose a maxChi2DoF corresponding to a significance level of at most 0.05
    # (note this is actually a significance not a chi2 number)
    maxChi2DoF = pexConfig.Field(
        dtype=float, default=0.05,
        doc="Maximum Chi2/DoF of fit to be considered a dipole")

    verbose = pexConfig.Field(
        dtype=bool, default=False,
        doc="be verbose; this is slow")


class DipoleFitTask(measBase.SingleFrameMeasurementTask):
    """This is a subclass of SingleFrameMeasurementTask which can accept
    three input images in its run() method. Because it subclasses
    SingleFrameMeasurementTask, and calls
    SingleFrameMeasurementTask.run() from its run() method, it still
    can be used identically to a standard SingleFrameMeasurementTask.

    """

    ConfigClass = DipoleFitConfig
    _DefaultName = "ip_diffim_DipoleFit"

    def __init__(self, schema, algMetadata=None, dpFitConfig=None, **kwds):

        measBase.SingleFrameMeasurementTask.__init__(self, schema, algMetadata, **kwds)

        self.dpFitConfig = dpFitConfig
        if self.dpFitConfig is None:
            self.dpFitConfig = DipoleFitConfig()
        self.dipoleFitter = DipoleFitPlugin(self.dpFitConfig, name=self._DefaultName,
                                            schema=schema, metadata=algMetadata)

    def run(self, sources, exposure, posImage=None, negImage=None, **kwds):
        """!Run dipole measurement and classification

        @param sources       diaSources that will be measured using dipole measurement
        @param exposure      Exposure on which the diaSources were detected
        @param posImage      "Positive" exposure from which the template was subtracted
        @param negImage      "Negative" exposure which was subtracted from the posImage
        @param **kwds        Sent to SingleFrameMeasurementTask
        """

        measBase.SingleFrameMeasurementTask.run(self, sources, exposure, **kwds)

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
    """Subclass of SingleFramePlugin which can accept three input images
    in its measure() method and fits dipoles to all merged (two-peak)
    footprints in a diffim/pos-im/neg-im simultaneously. The meat of
    the fitting routines are in the class DipoleFitAlgorithm.

    The motivation behind this plugin and the necessity for including more than
    one exposure are documented in DMTN-007 (http://dmtn-007.readthedocs.org).

    This class is named ip_diffim_DipoleFit so that it may be used alongside
    the existing ip_diffim_DipoleMeasurement classes until such a time as those
    are deemed to be replaceable by this.

    """

    ConfigClass = DipoleFitConfig

    FAILURE_EDGE = 1   # too close to the edge
    FAILURE_FIT = 2    # failure in the fitting
    FAILURE_NOT_DIPOLE = 4  # input source is not a putative dipole to begin with

    @classmethod
    def getExecutionOrder(cls):
        # algorithms that require both getShape() and getCentroid(),
        # in addition to a Footprint and its Peaks
        return cls.FLUX_ORDER

    @classmethod
    def getTransformClass(cls):
        return DipoleFitTransform

    def __init__(self, config, name, schema, metadata):
        measBase.SingleFramePlugin.__init__(self, config, name, schema, metadata)

        self.log = lsst.pex.logging.Log(lsst.pex.logging.Log.getDefaultLog(),
                                        'lsst.ip.diffim.DipoleFitPlugin', lsst.pex.logging.Log.INFO)

        self._setupSchema(config, name, schema, metadata)

    def _setupSchema(self, config, name, schema, metadata):
        # Get a FunctorKey that can quickly look up the "blessed" centroid value.
        self.centroidKey = afwTable.Point2DKey(schema["slot_Centroid"])

        # Add some fields for our outputs, and save their Keys.
        # Use setattr() to programmatically set the pos/neg named attributes to values, e.g.
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
        """Perform the non-linear least squares minimization. The main
        functionality of this routine was placed outside of this
        plugin (into `DipoleFitAlgorithm.fitDipole()`) so that
        `DipoleFitAlgorithm.fitDipole()` can be called separately for
        testing (see `tests/testDipoleFitter.py`)

        """
        pks = measRecord.getFootprint().getPeaks()
        if len(pks) <= 1:  # not a dipole for our analysis
            self.fail(measRecord, measBase.MeasurementError('not a dipole', self.FAILURE_NOT_DIPOLE))

        try:
            alg = DipoleFitAlgorithm(exposure, posImage=posImage, negImage=negImage)
            result = alg.fitDipole(
                measRecord, rel_weight=self.config.relWeight,
                tol=self.config.tolerance,
                centroidRangeInSigma=self.config.centroidRange,
                fitBgGradient=self.config.fitBgGradient,
                separateNegParams=self.config.fitSeparateNegParams,
                verbose=self.config.verbose, display=False)
        except lsst.pex.exceptions.LengthError:
            raise measBase.MeasurementError('edge failure', self.FAILURE_EDGE)
        except Exception:
            self.fail(measRecord, measBase.MeasurementError('dipole fit failure', self.FAILURE_FIT))

        self.log.log(self.log.DEBUG, "Dipole fit result: %s" % str(result))

        if result.psfFitPosFlux <= 1.:   # usually around 0.1 -- the minimum flux allowed -- i.e. bad fit.
            self.fail(measRecord, measBase.MeasurementError('dipole fit failure', self.FAILURE_FIT))

        # add chi2, coord/flux uncertainties (TBD), dipole classification

        # Add the relevant values to the measRecord
        measRecord[self.posFluxKey] = result.psfFitPosFlux
        measRecord[self.posFluxSigmaKey] = result.psfFitSignaltoNoise   # to be changed to actual sigma!
        measRecord[self.posCentroidKeyX] = result.psfFitPosCentroidX
        measRecord[self.posCentroidKeyY] = result.psfFitPosCentroidY

        measRecord[self.negFluxKey] = result.psfFitNegFlux
        measRecord[self.negFluxSigmaKey] = result.psfFitSignaltoNoise   # to be changed to actual sigma!
        measRecord[self.negCentroidKeyX] = result.psfFitNegCentroidX
        measRecord[self.negCentroidKeyY] = result.psfFitNegCentroidY

        # Dia source flux: average of pos+neg
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
        """Determine if source is classified as dipole (similar to
        orig. dipole classification task).
        """

        # First, does the total signal-to-noise surpass the minSn?
        passesSn = measRecord[self.signalToNoiseKey] > self.config.minSn

        # Second, are the pos/neg fluxes less than 1.0 and no more than 0.65 of the total flux?
        # By default this will never happen since posFlux = negFlux.
        passesFluxPos = (abs(measRecord[self.posFluxKey]) /
                         (measRecord[self.fluxKey]*2.)) < self.config.maxFluxRatio
        passesFluxPos &= (abs(measRecord[self.posFluxKey]) >= 1.0)
        passesFluxNeg = (abs(measRecord[self.negFluxKey]) /
                         (measRecord[self.fluxKey]*2.)) < self.config.maxFluxRatio
        passesFluxNeg &= (abs(measRecord[self.negFluxKey]) >= 1.0)
        allPass = (passesSn and passesFluxPos and passesFluxNeg)  # and passesChi2)

        # Third, is it a good fit (chi2dof < 1)?
        # Use scipy's chi2 cumulative distrib to estimate significance
        # This doesn't really work since I don't trust the values in the variance plane (which
        #   affects the least-sq weights, which affects the resulting chi2).
        # But I'm going to keep this here for future use.
        if False:
            from scipy.stats import chi2
            ndof = result.psfFitChi2 / measRecord[self.chi2dofKey]
            significance = chi2.cdf(result.psfFitChi2, ndof)
            passesChi2 = significance < self.config.maxChi2DoF
            allPass = allPass and passesChi2

        if allPass:  # Note cannot pass `allPass` into the `measRecord.set()` call below...?
            measRecord.set(self.classificationFlagKey, True)
        else:
            measRecord.set(self.classificationFlagKey, False)

    # TBD: need to catch more exceptions, set correct flags.
    def fail(self, measRecord, error=None):
        measRecord.set(self.flagKey, True)
        if error is not None:
            if error.getFlagBit() == self.FAILURE_EDGE:
                measRecord.set(self.edgeFlagKey, True)
            if error.getFlagBit() == self.FAILURE_FIT:
                measRecord.set(self.flagKey, True)
            if error.getFlagBit() == self.FAILURE_NOT_DIPOLE:
                measRecord.set(self.flagKey, True)


class DipoleFitAlgorithm(object):
    """Lightweight class containing methods for fitting dipoles in
    diffims, used by DipoleFitPlugin.  This code is documented in
    DMTN-007.

    Below is a (somewhat incomplete) list of improvements
    that would be worth investigating, given the time:

    1. Initial fast test whether a background gradient needs to be fit
    2. Initial fast estimate of background gradient(s) params --
    perhaps using numpy.lstsq
    3. better estimate for staring flux when there's a strong gradient
    4. evaluate necessity for separate parameters for pos- and neg-
    images
    5. only fit background OUTSIDE footprint and dipole params INSIDE
    footprint?
    6. requires a new package `lmfit` -- investiate others?
    (astropy/scipy/iminuit?)
    7. account for PSFs that vary across the exposures

    """

    # This is just a private version number to sync with the ipython notebooks that I have been
    # using for algorithm development.
    __private_version__ = '0.0.3'

    # Create a namedtuple to hold all of the relevant output from the lmfit results
    resultsOutput = namedtuple('resultsOutput',
                               ['psfFitPosCentroidX', 'psfFitPosCentroidY',
                                'psfFitNegCentroidX', 'psfFitNegCentroidY', 'psfFitPosFlux',
                                'psfFitNegFlux', 'psfFitPosFluxSigma', 'psfFitNegFluxSigma',
                                'psfFitCentroidX', 'psfFitCentroidY', 'psfFitOrientation',
                                'psfFitSignaltoNoise', 'psfFitChi2', 'psfFitRedChi2'])

    def __init__(self, diffim, posImage=None, negImage=None):
        """Algorithm to run dipole measurement on a diaSource
        @param diffim Exposure on which the diaSources were detected
        @param posImage "Positive" exposure from which the template was subtracted
        @param negImage "Negative" exposure which was subtracted from the posImage
        """

        self.diffim = diffim
        self.posImage = posImage
        self.negImage = negImage

    def genBgGradientModel(self, in_x, b=None, x1=0., y1=0., xy=None, x2=0., y2=0.):
        """Generate gradient model (2-d array) with up to 2nd-order polynomial
        @param in_x Provides the input x,y grid upon which to compute the gradient
        @param b Intercept (0th-order parameter) for gradient. If None, do nothing (for speed).
        @param x1 X-slope (1st-order parameter) for gradient.
        @param y1 Y-slope (1st-order parameter) for gradient.
        @param xy X*Y coefficient (2nd-order parameter) for gradient.
        If None, do not compute 2nd order polynomial.
        @param x2 X**2 coefficient (2nd-order parameter) for gradient.
        @param y2 Y**2 coefficient (2nd-order parameter) for gradient.

        @return None, or 2-d numpy.array of width/height matching
        input bbox, containing computed gradient values.

        """

        gradient = None
        if b is not None:  # Don't fit for other gradient parameters if the intercept is not allowed.
            y, x = in_x[0, :], in_x[1, :]
            gradient = np.full_like(x, b, dtype='float64')
            if x1 is not None:
                gradient += x1 * x
            if y1 is not None:
                gradient += y1 * y
            if xy is not None:
                gradient += xy * (x * y)
            if x2 is not None:
                gradient += x2 * (x * x)
            if y2 is not None:
                gradient += y2 * (y * y)
        return gradient

    def genStarModel(self, bbox, psf, xcen, ycen, flux):
        """Generate model (2-d array) of a 'star' (single PSF) centered at given coordinates

        @param bbox Bounding box marking pixel coordinates for generated model
        @param psf Psf model used to generate the 'star'
        @param xcen Desired x-centroid of the 'star'
        @param ycen Desired y-centroid of the 'star'
        @param flux Desired flux of the 'star'

        @return 2-d numpy.array of width/height matching input bbox,
        containing PSF with given centroid and flux

        """

        # Generate the psf image, normalize to flux
        psf_img = psf.computeImage(afwGeom.Point2D(xcen, ycen)).convertF()
        psf_img_sum = np.nansum(psf_img.getArray())
        psf_img *= (flux/psf_img_sum)

        # Clip the PSF image bounding box to fall within the footprint bounding box
        psf_box = psf_img.getBBox()
        psf_box.clip(bbox)
        psf_img = afwImage.ImageF(psf_img, psf_box, afwImage.PARENT)

        # Then actually crop the psf image.
        # Usually not necessary, but if the dipole is near the edge of the image...
        # Would be nice if we could compare original pos_box with clipped pos_box and
        #     see if it actually was clipped.
        p_Im = afwImage.ImageF(bbox)
        tmpSubim = afwImage.ImageF(p_Im, psf_box, afwImage.PARENT)
        tmpSubim += psf_img

        return p_Im

    @staticmethod
    def genDipoleModel(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=None,
                       b=None, x1=None, y1=None, xy=None, x2=None, y2=None,
                       bNeg=None, x1Neg=None, y1Neg=None, xyNeg=None, x2Neg=None, y2Neg=None,
                       debug=False, **kwargs):
        """Generate dipole model with given parameters.

        This is the functor whose sum-of-squared difference from data
        is minimized by `lmfit`. Thus it must be static. However, it
        just defers to `self.genDipoleModelImpl()`, where `self` comes
        out of kwargs['myself'].

        @see genDipoleModelImpl

        """
        algObject = kwargs.pop('algObject')
        return algObject.genDipoleModelImpl(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=fluxNeg,
                                            b=b, x1=x1, y1=y1, xy=xy, x2=x2, y2=y2,
                                            bNeg=bNeg, x1Neg=x1Neg, y1Neg=y1Neg, xyNeg=xyNeg,
                                            x2Neg=x2Neg, y2Neg=y2Neg, debug=debug, **kwargs)

    def genDipoleModelImpl(self, x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=None,
                       b=None, x1=None, y1=None, xy=None, x2=None, y2=None,
                       bNeg=None, x1Neg=None, y1Neg=None, xyNeg=None, x2Neg=None, y2Neg=None,
                       debug=False, **kwargs):

        """Generate dipole model with given parameters. This is the functor
        whose sum-of-squared difference from data is minimized by
        `lmfit`.

        @param x Input independent variable. Used here as the grid on
        which to compute the background gradient model.
        @param flux Desired flux of the positive lobe of the dipole
        @param xcenPos Desired x-centroid of the positive lobe of the dipole
        @param ycenPos Desired y-centroid of the positive lobe of the dipole
        @param xcenNeg Desired x-centroid of the negative lobe of the dipole
        @param ycenNeg Desired y-centroid of the negative lobe of the dipole
        @param fluxNeg Desired flux of the negative lobe of the dipole, set to 'flux' if None
        @param b, x1, y1, xy, x2, y2 Gradient parameters for positive lobe.
        @see `genBgGradientModel`.
        @param bNeg, x1Neg, y1Neg, xyNeg, x2Neg, y2Neg Gradient
        parameters for negative lobe.  @see `genBgGradientModel`. Set
        to the corresponding positive values if None.

        @param **kwargs Keyword arguments passed through `lmfit` and
        used by this functor. These must include:
           @param psf Psf model used to generate the 'star'
           @param rel_weight Used to signify if positive/negative images are to be
           included (!= 0. if yes)
           @param bbox Bounding box containing region to be modelled

        @return numpy.array of width/height matching input bbox,
        containing dipole model with given centroids and flux(es). If
        rel_weight = 0., this is a 2-d array with dimensions matching
        those of bbox; otherwise a stack of three such arrays,
        representing the dipole (diffim), positive and negative images
        respectively.

        """

        psf = kwargs.get('psf')
        rel_weight = kwargs.get('rel_weight')  # if > 0, we're including pre-sub. images
        fp = kwargs.get('footprint')
        bbox = fp.getBBox()

        if fluxNeg is None:
            fluxNeg = flux

        if debug:
            print('%.2f %.2f %.2f %.2f %.2f %.2f' % (flux, fluxNeg, xcenPos, ycenPos, xcenNeg, ycenNeg))
            if x1 is not None:
                print('     %.2f %.2f %.2f' % (b, x1, y1))
            if xy is not None:
                print('     %.2f %.2f %.2f' % (xy, x2, y2))

        posIm = self.genStarModel(bbox, psf, xcenPos, ycenPos, flux)
        negIm = self.genStarModel(bbox, psf, xcenNeg, ycenNeg, fluxNeg)

        in_x = x
        if in_x is None:  # use the footprint to generate the input grid
            y, x = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
            in_x = np.array([x, y]) * 1.

        gradient = self.genBgGradientModel(in_x, b, x1, y1, xy, x2, y2)
        gradientNeg = gradient
        if bNeg is not None:
            gradientNeg = self.genBgGradientModel(in_x, bNeg, x1Neg, y1Neg, xyNeg, x2Neg, y2Neg)

        if gradient is not None:
            posIm.getArray()[:, :] += gradient
            negIm.getArray()[:, :] += gradientNeg

        # Generate the diffIm model
        diffIm = afwImage.ImageF(bbox)
        diffIm += posIm
        diffIm -= negIm

        zout = diffIm.getArray()
        if rel_weight > 0.:
            zout = np.append([zout], [posIm.getArray(), negIm.getArray()], axis=0)

        return zout

    def fitDipoleImpl(self, source, tol=1e-7, rel_weight=0.5,
                      fitBgGradient=True, bgGradientOrder=1, centroidRangeInSigma=5.,
                      separateNegParams=True, verbose=False, display=False):
        """Fit a dipole model to an input difference image (actually,
        subimage bounded by the input source's footprint) and
        optionally constrain the fit using the pre-subtraction images
        posImage and negImage.

        @see `fitDipole()`
        @return `lmfit.MinimizerResult` object containing the fit
        parameters and other information.

        """

        # Only import lmfit if someone wants to use the new DipoleFitAlgorithm.
        import lmfit

        fp = source.getFootprint()
        bbox = fp.getBBox()
        subim = afwImage.MaskedImageF(self.diffim.getMaskedImage(), bbox, afwImage.PARENT)

        z = diArr = subim.getArrays()[0]
        weights = 1. / subim.getArrays()[2]  # get the weights (=1/variance)
        if self.posImage is not None and rel_weight > 0.:
            posSubim = afwImage.MaskedImageF(self.posImage.getMaskedImage(), bbox, afwImage.PARENT)
            negSubim = afwImage.MaskedImageF(self.negImage.getMaskedImage(), bbox, afwImage.PARENT)
            z = np.append([z], [posSubim.getArrays()[0],
                                negSubim.getArrays()[0]], axis=0)
            # Weight the pos/neg images by rel_weight relative to the diffim
            weights = np.append([weights], [1. / posSubim.getArrays()[2] * rel_weight,
                                            1. / negSubim.getArrays()[2] * rel_weight], axis=0)

        psfSigma = self.diffim.getPsf().computeShape().getDeterminantRadius()

        # Create the lmfit model (lmfit uses scipy 'leastsq' option by default - Levenberg-Marquardt)
        gmod = lmfit.Model(DipoleFitAlgorithm.genDipoleModel, verbose=verbose)

        # Add the constraints for centroids, fluxes.
        # starting constraint - near centroid of footprint
        fpCentroid = np.array([fp.getCentroid().getX(), fp.getCentroid().getY()])
        cenNeg = cenPos = fpCentroid

        pks = fp.getPeaks()
        if len(pks) >= 1:
            cenPos = pks[0].getF()    # if individual (merged) peaks were detected, use those
        if len(pks) >= 2:
            cenNeg = pks[1].getF()

        # For close/faint dipoles the starting locs (min/max) might be way off, let's help them a bit.
        # First assume dipole is not separated by more than 5*psfSigma.
        centroidRange = psfSigma * centroidRangeInSigma

        # Note - this may be a cheat to assume the dipole is centered in center of the footprint.
        if np.sum(np.sqrt((np.array(cenPos) - fpCentroid)**2.)) > centroidRange:
            cenPos = fpCentroid
        if np.sum(np.sqrt((np.array(cenNeg) - fpCentroid)**2.)) > centroidRange:
            cenPos = fpCentroid

        # parameter hints/constraints: https://lmfit.github.io/lmfit-py/model.html#model-param-hints-section
        # might make sense to not use bounds -- see http://lmfit.github.io/lmfit-py/bounds.html
        # also see this discussion -- https://github.com/scipy/scipy/issues/3129
        gmod.set_param_hint('xcenPos', value=cenPos[0],
                            min=cenPos[0]-centroidRange, max=cenPos[0]+centroidRange)
        gmod.set_param_hint('ycenPos', value=cenPos[1],
                            min=cenPos[1]-centroidRange, max=cenPos[1]+centroidRange)
        gmod.set_param_hint('xcenNeg', value=cenNeg[0],
                            min=cenNeg[0]-centroidRange, max=cenNeg[0]+centroidRange)
        gmod.set_param_hint('ycenNeg', value=cenNeg[1],
                            min=cenNeg[1]-centroidRange, max=cenNeg[1]+centroidRange)

        # Estimate starting flux. This strongly affects runtime performance so we want to make it close.
        # I tried many possibilities, the best one seems to be just using sum(abs(diffIm))/2.
        # I am leaving the other tries here (commented out) for posterity and possible future improvements...

        # Value to convert peak value to total flux based on flux within psf
        # psfImg = diffim.getPsf().computeImage()
        # pkToFlux = np.nansum(psfImg.getArray()) / diffim.getPsf().computePeak()

        # bg = np.nanmedian(diArr) # Compute the dipole background (probably very close to zero)
        # startingPk = np.nanmax(diArr) - bg   # use the dipole peak for an estimate.
        # posFlux, negFlux = startingPk * pkToFlux, -startingPk * pkToFlux

        # if len(pks) >= 1:
        #     posFlux = pks[0].getPeakValue() * pkToFlux
        # if len(pks) >= 2:
        #     negFlux = pks[1].getPeakValue() * pkToFlux

        # Use the (flux under the dipole)*5 for an estimate.
        # Lots of testing showed that having startingFlux be too high was better than too low.
        startingFlux = np.nansum(np.abs(diArr) - np.nanmedian(np.abs(diArr))) * 5.
        posFlux = negFlux = startingFlux

        # TBD: set max. flux limit?
        gmod.set_param_hint('flux', value=posFlux, min=0.1)

        if separateNegParams:
            # TBD: set max negative lobe flux limit?
            gmod.set_param_hint('fluxNeg', value=np.abs(negFlux), min=0.1)

        # Fixed parameters (dont fit for them if there are no pre-sub images or no gradient fit requested):
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

        y, x = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
        in_x = np.array([x, y]).astype(np.float)

        # I'm not sure about the variance planes in the diffim (or convolved pre-sub. images
        # for that matter) so for now, let's just do an un-weighted least-squares fit
        # (override weights computed above).
        weights = 1.
        if self.posImage is not None and rel_weight > 0.:
            weights = np.array([np.ones_like(diArr), np.ones_like(diArr)*rel_weight,
                                np.ones_like(diArr)*rel_weight])

        # Note that although we can, we're not required to set initial values for params here,
        # since we set their param_hint's above.
        # add "method" param to not use 'leastsq' (==levenberg-marquardt), e.g. "method='nelder'"
        result = gmod.fit(z, weights=weights, x=in_x,
                          verbose=verbose,
                          fit_kws={'ftol': tol, 'xtol': tol, 'gtol': tol, 'maxfev': 250},  # see scipy docs
                          psf=self.diffim.getPsf(),  # hereon: kwargs that get passed to genDipoleModel()
                          rel_weight=rel_weight,
                          footprint=fp,
                          algObject=self)

        # Probably never wanted - also this takes a long time (longer than the fit!)
        # This is how to get confidence intervals out:
        #    https://lmfit.github.io/lmfit-py/confidence.html and
        #    http://cars9.uchicago.edu/software/python/lmfit/model.html
        if verbose:  # the ci_report() seems to fail if neg params are constrained -- TBD why.
            print(result.fit_report(show_correl=False))
            if separateNegParams:
                print(result.ci_report())

        # Display images, model fits and residuals (currently uses matplotlib display functions)
        if display:
            self.displayFitResults(result, fp)

        return result

    def fitDipole(self, source, tol=1e-7, rel_weight=0.1,
                  fitBgGradient=True, centroidRangeInSigma=5., separateNegParams=True,
                  bgGradientOrder=1, verbose=False, display=False, return_fitObj=False):
        """Wrapper around `fitDipoleImpl()` which performs the fit of a dipole
        model to an input difference image (actually, subimage bounded
        by the input source's footprint) and optionally constrain the
        fit using the pre-subtraction images posImage and
        negImage. Wraps the output into a `resultsOutput` object after
        computing additional statistics such as orientation and SNR.

        @param source Record containing the (merged) dipole source footprint detected on the diffim
        @param tol Tolerance parameter for scipy.leastsq() optimization
        @param rel_weight Weighting of posImage/negImage relative to the diffim in the fit
        @param fitBgGradient Fit linear background gradient in posImage/negImage?
        @param bgGradientOrder Desired polynomial order of background gradient (allowed are [0,1,2])
        @param centroidRangeInSigma Allowed window of centroid parameters relative to peak in input source footprint
        @param separateNegParams Fit separate parameters to the flux and background gradient in the negative images?
        @param verbose Be verbose
        @param display Display input data, best fit model(s) and residuals in a matplotlib window.
        @param return_fitObj In addition to the `resultsOutput` object, also returns the
        `lmfit.MinimizerResult` object for debugging.

        @return resultsOutput object containing the fit parameters and other information.
        @return lmfit.MinimizerResult object if `return_fitObj` is True.

        """

        fitResult = self.fitDipoleImpl(
            source, tol=tol, rel_weight=rel_weight, fitBgGradient=fitBgGradient,
            centroidRangeInSigma=centroidRangeInSigma, separateNegParams=separateNegParams,
            bgGradientOrder=bgGradientOrder, verbose=verbose, display=display)

        fitParams = fitResult.best_values
        if fitParams['flux'] <= 1.:   # usually around 0.1 -- the minimun flux allowed -- i.e. bad fit.
            out = DipoleFitAlgorithm.resultsOutput(
                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                np.nan, np.nan, np.nan)
            if return_fitObj:  # for debugging
                return out, fitResult
            return out

        centroid = ((fitParams['xcenPos'] + fitParams['xcenNeg']) / 2.,
                    (fitParams['ycenPos'] + fitParams['ycenNeg']) / 2.)
        dx, dy = fitParams['xcenPos'] - fitParams['xcenNeg'], fitParams['ycenPos'] - fitParams['ycenNeg']
        angle = np.arctan2(dy, dx) / np.pi * 180.   # convert to degrees (should keep as rad?)

        # Exctract flux value, compute signalToNoise from flux/variance_within_footprint
        # Also extract the stderr of flux estimate.
        def computeSumVariance(exposure, footprint):
            box = footprint.getBBox()
            subim = afwImage.MaskedImageF(exposure.getMaskedImage(), box, afwImage.PARENT)
            return np.sqrt(np.nansum(subim.getArrays()[1][:, :]))

        fluxVal = fluxVar = fitParams['flux']
        fluxErr = fluxErrNeg = fitResult.params['flux'].stderr
        if self.posImage is not None:
            fluxVar = computeSumVariance(self.posImage, source.getFootprint())
        else:
            fluxVar = computeSumVariance(self.diffim, source.getFootprint())

        fluxValNeg, fluxVarNeg = fluxVal, fluxVar
        if separateNegParams:
            fluxValNeg = fitParams['fluxNeg']
            fluxErrNeg = fitResult.params['fluxNeg'].stderr
        if self.negImage is not None:
            fluxVarNeg = computeSumVariance(self.negImage, source.getFootprint())

        try:
            signalToNoise = np.sqrt((fluxVal/fluxVar)**2 + (fluxValNeg/fluxVarNeg)**2)
        except:  # catch divide by zero - should never happen.
            signalToNoise = np.nan

        out = DipoleFitAlgorithm.resultsOutput(
            fitParams['xcenPos'], fitParams['ycenPos'], fitParams['xcenNeg'], fitParams['ycenNeg'],
            fluxVal, -fluxValNeg, fluxErr, fluxErrNeg, centroid[0], centroid[1], angle,
            signalToNoise, fitResult.chisqr, fitResult.redchi)

        if return_fitObj:  # for debugging
            return out, fitResult
        return out

    def displayFitResults(self, result, footprint):
        """Display data, model fits and residuals (currently uses matplotlib display functions)."""

        try:
            plt = DipolePlotUtils.importMatplotlib()
            if not plt:
                return

            z = result.data
            fit = result.best_fit
            bbox = footprint.getBBox()
            extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
            if z.shape[0] == 3:
                fig = plt.figure(figsize=(8, 8))
                for i in range(3):
                    plt.subplot(3, 3, i*3+1)
                    DipolePlotUtils.display2dArray(z[i, :], 'Data', True, extent=extent)
                    plt.subplot(3, 3, i*3+2)
                    DipolePlotUtils.display2dArray(fit[i, :], 'Model', True, extent=extent)
                    plt.subplot(3, 3, i*3+3)
                    DipolePlotUtils.display2dArray(z[i, :] - fit[i, :], 'Residual', True, extent=extent)
                return fig
            else:
                fig = plt.figure(figsize=(8, 2.5))
                plt.subplot(1, 3, 1)
                DipolePlotUtils.display2dArray(z, 'Data', True, extent=extent)
                plt.subplot(1, 3, 2)
                DipolePlotUtils.display2dArray(fit, 'Model', True, extent=extent)
                plt.subplot(1, 3, 3)
                DipolePlotUtils.display2dArray(z - fit, 'Residual', True, extent=extent)
                return fig
        except Exception as err:
            print('Uh oh! need matplotlib to use these funcs', err)
            pass


######### UTILITIES FUNCTIONS -- TBD WHERE THEY ULTIMATELY END UP ####

class DipolePlotUtils():
    """Utility class containing methods for displaying dipoles/footprints
    in difference images, mostly used for debugging.

    """

    @classmethod
    def importMatplotlib(cls):
        """!Prefer to import matplotlib.pyplot only when needed, allow for not
        available.
        @return the imported module if available, else False.

        """
        try:
            import matplotlib.pyplot as plt
            return plt
        except ImportError as err:
            print('Unable to import matplotlib: %s' % err)
            return False

    @classmethod
    def display2dArray(cls, arr, title='Data', showBars=True, extent=None):
        """Use matplotlib.pyplot.imshow() to display a 2-D array.

        @param arr The 2-D array to display
        @param title Optional title to display
        @param showBars Show grey-scale bar alongside
        @param extent If not None, a 4-tuple giving the bounding box coordinates of the array

        @return matplotlib.pyplot.figure dispaying the image
        """
        plt = cls.importMatplotlib()
        if not plt:
            return

        fig = plt.imshow(arr, origin='lower', interpolation='none', cmap='gray',
                         extent=extent)
        plt.title(title)
        if showBars:
            plt.colorbar(fig, cmap='gray')
        return fig

    @classmethod
    def displayImage(cls, image, showBars=True, width=8, height=2.5):
        """Use matplotlib.pyplot.imshow() to display an afw.image.Image within
        its bounding box.

        @see display2dArray() for params not listed here.
        @param image The image to display
        @param width The width of the display (inches)
        @param height The height of the display (inches)

        @return matplotlib.pyplot.figure dispaying the image

        """
        plt = cls.importMatplotlib()
        if not plt:
            return

        fig = plt.figure(figsize=(width, height))
        bbox = image.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
        plt.subplot(1, 3, 1)
        ma = image.getArray()
        cls.display2dArray(ma, title='Data', showBars=showBars, extent=extent)
        return fig

    @classmethod
    def displayImages(cls, images, showBars=True, width=8, height=2.5):
        """Use matplotlib.pyplot.imshow() to display up to three
        afw.image.Images alongside each other, each within its
        bounding box.

        @see displayImage() for params not listed here.
        @param images tuple of up to three images to display

        @return matplotlib.pyplot.figure dispaying the image

        """
        plt = cls.importMatplotlib()
        if not plt:
            return

        fig = plt.figure(figsize=(width, height))
        for i, image in enumerate(images):
            bbox = image.getBBox()
            extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
            plt.subplot(1, len(images), i+1)
            ma = image.getArray()
            cls.display2dArray(ma, title='Data', showBars=showBars, extent=extent)
        return fig

    @classmethod
    def displayMaskedImage(cls, maskedImage, showMasks=True, showVariance=False, showBars=True, width=8,
                           height=2.5):
        """Use matplotlib.pyplot.imshow() to display a afwImage.MaskedImageF,
        alongside its masks and variance plane

        !see displayImage() for params not listed here

        @param maskedImage MaskedImage to display
        @param showMasks Display the MaskedImage's masks
        @param showVariance Display the MaskedImage's variance plane

        @return matplotlib.pyplot.figure dispaying the image

        """
        plt = cls.importMatplotlib()
        if not plt:
            return

        fig = plt.figure(figsize=(width, height))
        bbox = maskedImage.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
        plt.subplot(1, 3, 1)
        ma = maskedImage.getArrays()
        cls.display2dArray(ma[0], title='Data', showBars=showBars, extent=extent)
        if showMasks:
            plt.subplot(1, 3, 2)
            cls.display2dArray(ma[1], title='Masks', showBars=showBars, extent=extent)
        if showVariance:
            plt.subplot(1, 3, 3)
            cls.display2dArray(ma[2], title='Variance', showBars=showBars, extent=extent)
        return fig

    @classmethod
    def displayExposure(cls, exposure, showMasks=True, showVariance=False, showPsf=False, showBars=True,
                        width=8, height=2.5):
        """Use matplotlib.pyplot.imshow() to display a afw.image.Exposure,
        including its masks and variance plane and optionally its Psf.

        @see displayMaskedImage() for params not listed here

        @param exposure Exposure to display
        @param showPsf Display the exposure's Psf

        @return matplotlib.pyplot.figure dispaying the image

        """
        plt = cls.importMatplotlib()
        if not plt:
            return

        fig = cls.displayMaskedImage(exposure.getMaskedImage(), showMasks,
                                     showVariance=not showPsf,
                                     showBars=showBars, width=width, height=height)
        if showPsf:
            plt.subplot(1, 3, 3)
            psfIm = exposure.getPsf().computeImage()
            bbox = psfIm.getBBox()
            extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())
            cls.display2dArray(psfIm.getArray(), title='PSF', showBars=showBars, extent=extent)
        return fig

    @classmethod
    def displayCutouts(cls, source, exposure, posImage=None, negImage=None, asHeavyFootprint=False, title=''):
        """Use matplotlib.pyplot.imshow() to display cutouts within up to
        three afw.image.Exposure's, given by an input SourceRecord.

        @param source afw.table.SourceRecord defining the footprint to extract and display
        @param exposure afw.image.Exposure from which to extract the cutout to display
        @param posImage Second exposure from which to extract the cutout to display
        @param negImage Third exposure from which to extract the cutout to display
        @param asHeavyFootprint Display the cutouts as
        afw.detection.HeavyFootprint, with regions outside the
        footprint removed

        @return matplotlib.pyplot.figure dispaying the image

        """
        plt = cls.importMatplotlib()
        if not plt:
            return

        fp = source.getFootprint()
        bbox = fp.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())

        fig = plt.figure(figsize=(8, 2.5))
        if not asHeavyFootprint:
            subexp = afwImage.ImageF(exposure.getMaskedImage().getImage(), bbox, afwImage.PARENT)
        else:
            hfp = afwDet.HeavyFootprintF(fp, exposure.getMaskedImage())
            subexp = cls.getHeavyFootprintSubimage(hfp)
        plt.subplot(1, 3, 1)
        cls.display2dArray(subexp.getArray(), title=title+' Diffim', extent=extent)
        if posImage is not None:
            if not asHeavyFootprint:
                subexp = afwImage.ImageF(posImage.getMaskedImage().getImage(), bbox, afwImage.PARENT)
            else:
                hfp = afwDet.HeavyFootprintF(fp, posImage.getMaskedImage())
                subexp = cls.getHeavyFootprintSubimage(hfp)
            plt.subplot(1, 3, 2)
            cls.display2dArray(subexp.getArray(), title=title+' Pos', extent=extent)
        if negImage is not None:
            if not asHeavyFootprint:
                subexp = afwImage.ImageF(negImage.getMaskedImage().getImage(), bbox, afwImage.PARENT)
            else:
                hfp = afwDet.HeavyFootprintF(fp, negImage.getMaskedImage())
                subexp = cls.getHeavyFootprintSubimage(hfp)
            plt.subplot(1, 3, 3)
            cls.display2dArray(subexp.getArray(), title=title+' Neg', extent=extent)
        return fig

    @classmethod
    def makeHeavyCatalog(cls, catalog, exposure, verbose=False):
        """Turn all footprints in a catalog into heavy footprints."""
        for i, source in enumerate(catalog):
            fp = source.getFootprint()
            if not fp.isHeavy():
                if verbose:
                    print(i, 'not heavy => heavy')
                hfp = afwDet.HeavyFootprintF(fp, exposure.getMaskedImage())
                source.setFootprint(hfp)

        return catalog

    @classmethod
    def getHeavyFootprintSubimage(fp, badfill=np.nan):
        """Extract the image from a heavyFootprint as an ImageF."""
        hfp = afwDet.HeavyFootprintF_cast(fp)
        bbox = hfp.getBBox()

        subim2 = afwImage.ImageF(bbox, badfill)  # set the mask to NA (can use 0. if desired)
        afwDet.expandArray(hfp, hfp.getImageArray(), subim2.getArray(), bbox.getCorners()[0])
        return subim2

    @classmethod
    def searchCatalog(catalog, x, y):
        """Search a catalog for a source whose footprint contains the given x, y pixel coordinates."""
        for i, s in enumerate(catalog):
            bbox = s.getFootprint().getBBox()
            if bbox.contains(afwGeom.Point2D(x, y)):
                print(i)
                return s
        return None

