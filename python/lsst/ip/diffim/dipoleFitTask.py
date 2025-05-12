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

import math
import logging
import numpy as np
import warnings

import lsst.afw.image as afwImage
import lsst.meas.base as measBase
import lsst.afw.detection as afwDet
import lsst.geom as geom
import lsst.pex.exceptions as pexExcept
import lsst.pex.config as pexConfig
from lsst.pipe.base import Struct
from lsst.utils.timer import timeMethod

__all__ = ("DipoleFitTask", "DipoleFitPlugin", "DipoleFitTaskConfig", "DipoleFitPluginConfig",
           "DipoleFitAlgorithm")


# Create a new measurement task (`DipoleFitTask`) that can handle all other SFM tasks but can
# pass a separate pos- and neg- exposure/image to the `DipoleFitPlugin`s `run()` method.


class DipoleFitPluginConfig(measBase.SingleFramePluginConfig):
    """Configuration for DipoleFitPlugin
    """

    fitAllDiaSources = pexConfig.Field(
        dtype=bool, default=False,
        doc="""Attempte dipole fit of all diaSources (otherwise just the ones consisting of overlapping
        positive and negative footprints)""")

    maxSeparation = pexConfig.Field(
        dtype=float, default=5.,
        doc="Assume dipole is not separated by more than maxSeparation * psfSigma")

    relWeight = pexConfig.Field(
        dtype=float, default=0.5,
        doc="""Relative weighting of pre-subtraction images (higher -> greater influence of pre-sub.
        images on fit)""")

    tolerance = pexConfig.Field(
        dtype=float, default=1e-7,
        doc="Fit tolerance")

    fitBackground = pexConfig.Field(
        dtype=int, default=1,
        doc="Set whether and how to fit for linear gradient in pre-sub. images. Possible values:"
            "0: do not fit background at all"
            "1 (default): pre-fit the background using linear least squares and then do not fit it as part"
            "of the dipole fitting optimization"
            "2: pre-fit the background using linear least squares (as in 1), and use the parameter"
            "estimates from that fit as starting parameters for an integrated re-fit of the background")

    fitSeparateNegParams = pexConfig.Field(
        dtype=bool, default=False,
        doc="Include parameters to fit for negative values (flux, gradient) separately from pos.")

    # Config params for classification of detected diaSources as dipole or not
    minSn = pexConfig.Field(
        dtype=float, default=math.sqrt(2) * 5.0,
        doc="Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole")

    maxFluxRatio = pexConfig.Field(
        dtype=float, default=0.65,
        doc="Maximum flux ratio in either lobe to be considered a dipole")

    maxChi2DoF = pexConfig.Field(
        dtype=float, default=0.05,
        doc="""Maximum Chi2/DoF significance of fit to be considered a dipole.
        Default value means \"Choose a chi2DoF corresponding to a significance level of at most 0.05\"
        (note this is actually a significance, not a chi2 value).""")

    maxFootprintArea = pexConfig.Field(
        dtype=int, default=1_200,
        doc=("Maximum area for footprints before they are ignored as large; "
             "non-positive means no threshold applied"
             "Threshold chosen for HSC and DECam data, see DM-38741 for details."))


class DipoleFitTaskConfig(measBase.SingleFrameMeasurementConfig):

    def setDefaults(self):
        measBase.SingleFrameMeasurementConfig.setDefaults(self)

        self.plugins.names = ["base_SdssCentroid",
                              "ip_diffim_DipoleFit",
                              "base_CircularApertureFlux",
                              "base_PixelFlags",
                              "base_SkyCoord",
                              "base_PsfFlux",
                              ]
        # Only measure the apertures we need to report in the alert stream.
        self.plugins["base_CircularApertureFlux"].radii = [12.0]

        self.slots.calibFlux = None
        self.slots.modelFlux = None
        self.slots.gaussianFlux = None
        self.slots.shape = None
        # This will be switched to "ip_diffim_DipoleFit" as this task runs.
        self.slots.centroid = "base_SdssCentroid"
        self.doReplaceWithNoise = False


class DipoleFitTask(measBase.SingleFrameMeasurementTask):
    """A task that fits a dipole to a difference image, with an optional
    separate detection image.

    Because it subclasses SingleFrameMeasurementTask, and calls
    SingleFrameMeasurementTask.run() from its run() method, it still
    can be used identically to a standard SingleFrameMeasurementTask.
    """

    ConfigClass = DipoleFitTaskConfig
    _DefaultName = "dipoleFit"

    def __init__(self, schema, algMetadata=None, **kwargs):
        super().__init__(schema, algMetadata, **kwargs)

        # Enforce a specific plugin order, so that DipoleFit can fall back on
        # SdssCentroid for non-dipoles
        self.plugins_pre = self.plugins.copy()
        self.plugins_post = self.plugins.copy()
        self.plugins_pre.clear()
        self.plugins_pre["base_SdssCentroid"] = self.plugins["base_SdssCentroid"]
        self.plugins_post.pop("base_SdssCentroid")
        self.dipoleFit = self.plugins_post.pop("ip_diffim_DipoleFit")
        del self.plugins

    @timeMethod
    def run(self, sources, exposure, posExp=None, negExp=None, **kwargs):
        """Run dipole measurement and classification.

        Run SdssCentroid first, then switch the centroid slot, then DipoleFit
        then the rest; DipoleFit will fall back on SdssCentroid for sources
        not containing positive+negative peaks.

        Parameters
        ----------
        sources : `lsst.afw.table.SourceCatalog`
            ``diaSources`` that will be measured using dipole measurement.
        exposure : `lsst.afw.image.Exposure`
            The difference exposure on which the ``sources`` were detected.
            If neither ``posExp`` nor ``negExp`` are set, then the dipole is also
            fitted directly to this difference image.
        posExp : `lsst.afw.image.Exposure`, optional
            "Positive" exposure, typically a science exposure, or None if unavailable
            When `posExp` is `None`, will compute `posImage = exposure + negExp`.
        negExp : `lsst.afw.image.Exposure`, optional
            "Negative" exposure, typically a template exposure, or None if unavailable
            When `negExp` is `None`, will compute `negImage = posExp - exposure`.
        **kwargs
            Additional keyword arguments for `lsst.meas.base.sfm.SingleFrameMeasurementTask`.
        """
        # Run plugins in a very specific order, so DipoleFitPlugin has a
        # centroid to fall back on.
        self.plugins = self.plugins_pre
        super().run(sources, exposure, **kwargs)

        for source in sources:
            self.dipoleFit.measureDipoles(source, exposure, posExp, negExp)
        # Use the new DipoleFit outputs for subsequent measurements, now that
        # non-dipoles have been filled in with the earlier centroid values.
        sources.schema.getAliasMap().set("slot_Centroid", "ip_diffim_DipoleFit")

        self.plugins = self.plugins_post
        super().run(sources, exposure, **kwargs)


class DipoleModel:
    """Lightweight class containing methods for generating a dipole model for fitting
    to sources in diffims, used by DipoleFitAlgorithm.

    See also:
    `DMTN-007: Dipole characterization for image differencing  <https://dmtn-007.lsst.io>`_.
    """

    def __init__(self):
        import lsstDebug
        self.debug = lsstDebug.Info(__name__).debug
        self.log = logging.getLogger(__name__)

    def makeBackgroundModel(self, in_x, pars=None):
        """Generate gradient model (2-d array) with up to 2nd-order polynomial

        Parameters
        ----------
        in_x : `numpy.array`
            (2, w, h)-dimensional `numpy.array`, containing the
            input x,y meshgrid providing the coordinates upon which to
            compute the gradient. This will typically be generated via
            `_generateXYGrid()`. `w` and `h` correspond to the width and
            height of the desired grid.
        pars : `list` of `float`, optional
            Up to 6 floats for up
            to 6 2nd-order 2-d polynomial gradient parameters, in the
            following order: (intercept, x, y, xy, x**2, y**2). If `pars`
            is emtpy or `None`, do nothing and return `None` (for speed).

        Returns
        -------
        result : `None` or `numpy.array`
            return None, or 2-d numpy.array of width/height matching
            input bbox, containing computed gradient values.
        """

        # Don't fit for other gradient parameters if the intercept is not included.
        if (pars is None) or (len(pars) <= 0) or (pars[0] is None):
            return

        y, x = in_x[0, :], in_x[1, :]
        gradient = np.full_like(x, pars[0], dtype='float64')
        if len(pars) > 1 and pars[1] is not None:
            gradient += pars[1] * x
        if len(pars) > 2 and pars[2] is not None:
            gradient += pars[2] * y
        if len(pars) > 3 and pars[3] is not None:
            gradient += pars[3] * (x * y)
        if len(pars) > 4 and pars[4] is not None:
            gradient += pars[4] * (x * x)
        if len(pars) > 5 and pars[5] is not None:
            gradient += pars[5] * (y * y)

        return gradient

    def _generateXYGrid(self, bbox):
        """Generate a meshgrid covering the x,y coordinates bounded by bbox

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            input Bounding Box defining the coordinate limits

        Returns
        -------
        in_x : `numpy.array`
            (2, w, h)-dimensional numpy array containing the grid indexing over x- and
            y- coordinates
        """

        x, y = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
        in_x = np.array([y, x]).astype(np.float64)
        in_x[0, :] -= np.mean(in_x[0, :])
        in_x[1, :] -= np.mean(in_x[1, :])
        return in_x

    def _getHeavyFootprintSubimage(self, fp, badfill=np.nan, grow=0):
        """Extract the image from a ``~lsst.afw.detection.HeavyFootprint``
        as an `lsst.afw.image.ImageF`.

        Parameters
        ----------
        fp : `lsst.afw.detection.HeavyFootprint`
            HeavyFootprint to use to generate the subimage
        badfill : `float`, optional
            Value to fill in pixels in extracted image that are outside the footprint
        grow : `int`
            Optionally grow the footprint by this amount before extraction

        Returns
        -------
        subim2 : `lsst.afw.image.ImageF`
            An `~lsst.afw.image.ImageF` containing the subimage.
        """
        bbox = fp.getBBox()
        if grow > 0:
            bbox.grow(grow)

        subim2 = afwImage.ImageF(bbox, badfill)
        fp.getSpans().unflatten(subim2.array, fp.getImageArray(), bbox.getCorners()[0])
        return subim2

    def fitFootprintBackground(self, source, posImage, order=1):
        """Fit a linear (polynomial) model of given order (max 2) to the background of a footprint.

        Only fit the pixels OUTSIDE of the footprint, but within its bounding box.

        Parameters
        ----------
        source : `lsst.afw.table.SourceRecord`
            SourceRecord, the footprint of which is to be fit
        posImage : `lsst.afw.image.Exposure`
            The exposure from which to extract the footprint subimage
        order : `int`
            Polynomial order of background gradient to fit.

        Returns
        -------
        pars : `tuple` of `float`
            `tuple` of length (1 if order==0; 3 if order==1; 6 if order == 2),
            containing the resulting fit parameters
        """

        # TODO look into whether to use afwMath background methods -- see
        # http://lsst-web.ncsa.illinois.edu/doxygen/x_masterDoxyDoc/_background_example.html
        fp = source.getFootprint()
        bbox = fp.getBBox()
        bbox.grow(3)
        posImg = afwImage.ImageF(posImage.image, bbox, afwImage.PARENT)

        # This code constructs the footprint image so that we can identify the pixels that are
        # outside the footprint (but within the bounding box). These are the pixels used for
        # fitting the background.
        posHfp = afwDet.HeavyFootprintF(fp, posImage.getMaskedImage())
        posFpImg = self._getHeavyFootprintSubimage(posHfp, grow=3)

        isBg = np.isnan(posFpImg.array).ravel()

        data = posImg.array.ravel()
        data = data[isBg]
        B = data

        x, y = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
        x = x.astype(np.float64).ravel()
        x -= np.mean(x)
        x = x[isBg]
        y = y.astype(np.float64).ravel()
        y -= np.mean(y)
        y = y[isBg]
        b = np.ones_like(x, dtype=np.float64)

        M = np.vstack([b]).T  # order = 0
        if order == 1:
            M = np.vstack([b, x, y]).T
        elif order == 2:
            M = np.vstack([b, x, y, x**2., y**2., x*y]).T

        pars = np.linalg.lstsq(M, B, rcond=-1)[0]
        return pars

    def makeStarModel(self, bbox, psf, xcen, ycen, flux):
        """Generate a 2D image model of a single PDF centered at the given coordinates.

        Parameters
        ----------
        bbox : `lsst.geom.Box`
            Bounding box marking pixel coordinates for generated model
        psf : TODO: DM-17458
            Psf model used to generate the 'star'
        xcen : `float`
            Desired x-centroid of the 'star'
        ycen : `float`
            Desired y-centroid of the 'star'
        flux : `float`
            Desired flux of the 'star'

        Returns
        -------
        p_Im : `lsst.afw.image.Image`
            2-d stellar image of width/height matching input ``bbox``,
            containing PSF with given centroid and flux
        """

        # Generate the psf image, normalize to flux
        psf_img = psf.computeImage(geom.Point2D(xcen, ycen)).convertF()
        psf_img_sum = np.nansum(psf_img.array)
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

    def makeModel(self, x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=None,
                  b=None, x1=None, y1=None, xy=None, x2=None, y2=None,
                  bNeg=None, x1Neg=None, y1Neg=None, xyNeg=None, x2Neg=None, y2Neg=None,
                  **kwargs):
        """Generate dipole model with given parameters.

        This is the function whose sum-of-squared difference from data
        is minimized by `lmfit`.

        x : TODO: DM-17458
            Input independent variable. Used here as the grid on
            which to compute the background gradient model.
        flux : `float`
            Desired flux of the positive lobe of the dipole
        xcenPos, ycenPos : `float`
            Desired x,y-centroid of the positive lobe of the dipole
        xcenNeg, ycenNeg : `float`
            Desired x,y-centroid of the negative lobe of the dipole
        fluxNeg : `float`, optional
            Desired flux of the negative lobe of the dipole, set to 'flux' if None
        b, x1, y1, xy, x2, y2 : `float`
            Gradient parameters for positive lobe.
        bNeg, x1Neg, y1Neg, xyNeg, x2Neg, y2Neg : `float`, optional
            Gradient parameters for negative lobe.
            They are set to the corresponding positive values if None.

        **kwargs : `dict` [`str`]
            Keyword arguments passed through ``lmfit`` and
            used by this function. These must include:

            - ``psf`` Psf model used to generate the 'star'
            - ``rel_weight`` Used to signify least-squares weighting of posImage/negImage
                relative to diffim. If ``rel_weight == 0`` then posImage/negImage are ignored.
            - ``bbox`` Bounding box containing region to be modelled

        Returns
        -------
        zout : `numpy.array`
            Has width and height matching the input bbox, and
            contains the dipole model with given centroids and flux(es). If
            ``rel_weight`` = 0, this is a 2-d array with dimensions matching
            those of bbox; otherwise a stack of three such arrays,
            representing the dipole (diffim), positive, and negative images
            respectively.
        """

        psf = kwargs.get('psf')
        rel_weight = kwargs.get('rel_weight')  # if > 0, we're including pre-sub. images
        fp = kwargs.get('footprint')
        bbox = fp.getBBox()

        if fluxNeg is None:
            fluxNeg = flux

        self.log.debug('flux: %.2f fluxNeg: %.2f x+: %.2f x-: %.2f y+: %.2f y-: %.2f ',
                       flux, fluxNeg, xcenPos, xcenNeg, ycenPos, ycenNeg)
        if x1 is not None:
            self.log.debug('     b: %.2f x1: %.2f y1: %.2f', b, x1, y1)
        if xy is not None:
            self.log.debug('     xy: %.2f x2: %.2f y2: %.2f', xy, x2, y2)

        posIm = self.makeStarModel(bbox, psf, xcenPos, ycenPos, flux)
        negIm = self.makeStarModel(bbox, psf, xcenNeg, ycenNeg, fluxNeg)

        in_x = x
        if in_x is None:  # use the footprint to generate the input grid
            y, x = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
            in_x = np.array([x, y]) * 1.
            in_x[0, :] -= in_x[0, :].mean()  # center it!
            in_x[1, :] -= in_x[1, :].mean()

        if b is not None:
            gradient = self.makeBackgroundModel(in_x, (b, x1, y1, xy, x2, y2))

            # If bNeg is None, then don't fit the negative background separately
            if bNeg is not None:
                gradientNeg = self.makeBackgroundModel(in_x, (bNeg, x1Neg, y1Neg, xyNeg, x2Neg, y2Neg))
            else:
                gradientNeg = gradient

            posIm.array[:, :] += gradient
            negIm.array[:, :] += gradientNeg

        # Generate the diffIm model
        diffIm = afwImage.ImageF(bbox)
        diffIm += posIm
        diffIm -= negIm

        zout = diffIm.array
        if rel_weight > 0.:
            zout = np.append([zout], [posIm.array, negIm.array], axis=0)

        return zout


class DipoleFitAlgorithm:
    """Fit a dipole model using an image difference.

    See also:
    `DMTN-007: Dipole characterization for image differencing  <https://dmtn-007.lsst.io>`_.
    """

    # This is just a private version number to sync with the ipython notebooks that I have been
    # using for algorithm development.
    _private_version_ = '0.0.5'

    # Below is a (somewhat incomplete) list of improvements
    # that would be worth investigating, given the time:

    # todo 1. evaluate necessity for separate parameters for pos- and neg- images
    # todo 2. only fit background OUTSIDE footprint (DONE) and dipole params INSIDE footprint (NOT DONE)?
    # todo 3. correct normalization of least-squares weights based on variance planes
    # todo 4. account for PSFs that vary across the exposures (should be happening by default?)
    # todo 5. correctly account for NA/masks  (i.e., ignore!)
    # todo 6. better exception handling in the plugin
    # todo 7. better classification of dipoles (e.g. by comparing chi2 fit vs. monopole?)
    # todo 8. (DONE) Initial fast estimate of background gradient(s) params -- perhaps using numpy.lstsq
    # todo 9. (NOT NEEDED - see (2)) Initial fast test whether a background gradient needs to be fit
    # todo 10. (DONE) better initial estimate for flux when there's a strong gradient
    # todo 11. (DONE) requires a new package `lmfit` -- investiate others? (astropy/scipy/iminuit?)

    def __init__(self, diffim, posImage=None, negImage=None):
        """Algorithm to run dipole measurement on a diaSource

        Parameters
        ----------
        diffim : `lsst.afw.image.Exposure`
            Exposure on which the diaSources were detected
        posImage : `lsst.afw.image.Exposure`
            "Positive" exposure from which the template was subtracted
        negImage : `lsst.afw.image.Exposure`
            "Negative" exposure which was subtracted from the posImage
        """

        self.diffim = diffim
        self.posImage = posImage
        self.negImage = negImage
        self.psfSigma = None
        if diffim is not None:
            diffimPsf = diffim.getPsf()
            diffimAvgPos = diffimPsf.getAveragePosition()
            self.psfSigma = diffimPsf.computeShape(diffimAvgPos).getDeterminantRadius()

        self.log = logging.getLogger(__name__)

        import lsstDebug
        self.debug = lsstDebug.Info(__name__).debug

    def fitDipoleImpl(self, source, tol=1e-7, rel_weight=0.5,
                      fitBackground=1, bgGradientOrder=1, maxSepInSigma=5.,
                      separateNegParams=True, verbose=False):
        """Fit a dipole model to an input difference image.

        Actually, fits the subimage bounded by the input source's
        footprint) and optionally constrain the fit using the
        pre-subtraction images posImage and negImage.

        Parameters
        ----------
        source : TODO: DM-17458
            TODO: DM-17458
        tol : float, optional
            TODO: DM-17458
        rel_weight : `float`, optional
            TODO: DM-17458
        fitBackground : `int`, optional
            TODO: DM-17458
        bgGradientOrder : `int`, optional
            TODO: DM-17458
        maxSepInSigma : `float`, optional
            TODO: DM-17458
        separateNegParams : `bool`, optional
            TODO: DM-17458
        verbose : `bool`, optional
            TODO: DM-17458

        Returns
        -------
        result : `lmfit.MinimizerResult`
            return `lmfit.MinimizerResult` object containing the fit
            parameters and other information.
        """

        # Only import lmfit if someone wants to use the new DipoleFitAlgorithm.
        import lmfit

        fp = source.getFootprint()
        bbox = fp.getBBox()
        subim = afwImage.MaskedImageF(self.diffim.getMaskedImage(), bbox=bbox, origin=afwImage.PARENT)

        z = diArr = subim.image.array
        # Make sure we don't overwrite buffers.
        z = z.copy()
        weights = 1. / subim.variance.array  # get the weights (=1/variance)

        if rel_weight > 0. and ((self.posImage is not None) or (self.negImage is not None)):
            if self.negImage is not None:
                negSubim = afwImage.MaskedImageF(self.negImage.getMaskedImage(), bbox, origin=afwImage.PARENT)
            if self.posImage is not None:
                posSubim = afwImage.MaskedImageF(self.posImage.getMaskedImage(), bbox, origin=afwImage.PARENT)
            if self.posImage is None:  # no science image provided; generate it from diffim + negImage
                posSubim = subim.clone()
                posSubim += negSubim
            if self.negImage is None:  # no template provided; generate it from the posImage - diffim
                negSubim = posSubim.clone()
                negSubim -= subim

            z = np.append([z], [posSubim.image.array,
                                negSubim.image.array], axis=0)
            # Weight the pos/neg images by rel_weight relative to the diffim
            weights = np.append([weights], [1. / posSubim.variance.array * rel_weight,
                                            1. / negSubim.variance.array * rel_weight], axis=0)
        else:
            rel_weight = 0.  # a short-cut for "don't include the pre-subtraction data"

        # It seems that `lmfit` requires a static functor as its optimized method, which eliminates
        # the ability to pass a bound method or other class method. Here we write a wrapper which
        # makes this possible.
        def dipoleModelFunctor(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=None,
                               b=None, x1=None, y1=None, xy=None, x2=None, y2=None,
                               bNeg=None, x1Neg=None, y1Neg=None, xyNeg=None, x2Neg=None, y2Neg=None,
                               **kwargs):
            """Generate dipole model with given parameters.

            It simply defers to `modelObj.makeModel()`, where `modelObj` comes
            out of `kwargs['modelObj']`.
            """
            modelObj = kwargs.pop('modelObj')
            return modelObj.makeModel(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=fluxNeg,
                                      b=b, x1=x1, y1=y1, xy=xy, x2=x2, y2=y2,
                                      bNeg=bNeg, x1Neg=x1Neg, y1Neg=y1Neg, xyNeg=xyNeg,
                                      x2Neg=x2Neg, y2Neg=y2Neg, **kwargs)

        dipoleModel = DipoleModel()

        modelFunctor = dipoleModelFunctor  # dipoleModel.makeModel does not work for now.
        # Create the lmfit model (lmfit uses scipy 'leastsq' option by default - Levenberg-Marquardt)
        # We have to (later) filter out the nans by hand in our input to gmod.fit().
        # The only independent variable in the model is "x"; lmfit tries to
        # introspect variables and parameters from the function signature, but
        # gets it wrong for the model signature above.
        gmod = lmfit.Model(modelFunctor, independent_vars=["x"], verbose=verbose)

        # Add the constraints for centroids, fluxes.
        # starting constraint - near centroid of footprint
        fpCentroid = np.array([fp.getCentroid().getX(), fp.getCentroid().getY()])
        cenNeg = cenPos = fpCentroid

        pks = fp.getPeaks()

        if len(pks) >= 1:
            cenPos = pks[0].getF()    # if individual (merged) peaks were detected, use those
        if len(pks) >= 2:    # peaks are already sorted by centroid flux so take the most negative one
            cenNeg = pks[-1].getF()

        # For close/faint dipoles the starting locs (min/max) might be way off, let's help them a bit.
        # First assume dipole is not separated by more than 5*psfSigma.
        maxSep = self.psfSigma * maxSepInSigma

        # As an initial guess -- assume the dipole is close to the center of the footprint.
        if np.sum(np.sqrt((np.array(cenPos) - fpCentroid)**2.)) > maxSep:
            cenPos = fpCentroid
        if np.sum(np.sqrt((np.array(cenNeg) - fpCentroid)**2.)) > maxSep:
            cenPos = fpCentroid

        # parameter hints/constraints: https://lmfit.github.io/lmfit-py/model.html#model-param-hints-section
        # might make sense to not use bounds -- see http://lmfit.github.io/lmfit-py/bounds.html
        # also see this discussion -- https://github.com/scipy/scipy/issues/3129
        gmod.set_param_hint('xcenPos', value=cenPos[0],
                            min=cenPos[0]-maxSep, max=cenPos[0]+maxSep)
        gmod.set_param_hint('ycenPos', value=cenPos[1],
                            min=cenPos[1]-maxSep, max=cenPos[1]+maxSep)
        gmod.set_param_hint('xcenNeg', value=cenNeg[0],
                            min=cenNeg[0]-maxSep, max=cenNeg[0]+maxSep)
        gmod.set_param_hint('ycenNeg', value=cenNeg[1],
                            min=cenNeg[1]-maxSep, max=cenNeg[1]+maxSep)

        # Use the (flux under the dipole)*5 for an estimate.
        # Lots of testing showed that having startingFlux be too high was better than too low.
        startingFlux = np.nansum(np.abs(diArr) - np.nanmedian(np.abs(diArr))) * 5.
        posFlux = negFlux = startingFlux

        # TBD: set max. flux limit?
        gmod.set_param_hint('flux', value=posFlux, min=0.1)

        if separateNegParams:
            # TBD: set max negative lobe flux limit?
            gmod.set_param_hint('fluxNeg', value=np.abs(negFlux), min=0.1)

        # Fixed parameters (don't fit for them if there are no pre-sub images or no gradient fit requested):
        # Right now (fitBackground == 1), we fit a linear model to the background and then subtract
        # it from the data and then don't fit the background again (this is faster).
        # A slower alternative (fitBackground == 2) is to use the estimated background parameters as
        # starting points in the integrated model fit. That is currently not performed by default,
        # but might be desirable in some cases.
        bgParsPos = bgParsNeg = (0., 0., 0.)
        if ((rel_weight > 0.) and (fitBackground != 0) and (bgGradientOrder >= 0)):
            pbg = 0.
            bgFitImage = self.posImage if self.posImage is not None else self.negImage
            # Fit the gradient to the background (linear model)
            bgParsPos = bgParsNeg = dipoleModel.fitFootprintBackground(source, bgFitImage,
                                                                       order=bgGradientOrder)

            # Generate the gradient and subtract it from the pre-subtraction image data
            if fitBackground == 1:
                in_x = dipoleModel._generateXYGrid(bbox)
                pbg = dipoleModel.makeBackgroundModel(in_x, tuple(bgParsPos))
                z[1, :] -= pbg
                z[1, :] -= np.nanmedian(z[1, :])
                posFlux = np.nansum(z[1, :])
                gmod.set_param_hint('flux', value=posFlux*1.5, min=0.1)

                if separateNegParams and self.negImage is not None:
                    bgParsNeg = dipoleModel.fitFootprintBackground(source, self.negImage,
                                                                   order=bgGradientOrder)
                    pbg = dipoleModel.makeBackgroundModel(in_x, tuple(bgParsNeg))
                z[2, :] -= pbg
                z[2, :] -= np.nanmedian(z[2, :])
                if separateNegParams:
                    negFlux = np.nansum(z[2, :])
                    gmod.set_param_hint('fluxNeg', value=negFlux*1.5, min=0.1)

            # Do not subtract the background from the images but include the background parameters in the fit
            if fitBackground == 2:
                if bgGradientOrder >= 0:
                    gmod.set_param_hint('b', value=bgParsPos[0])
                    if separateNegParams:
                        gmod.set_param_hint('bNeg', value=bgParsNeg[0])
                if bgGradientOrder >= 1:
                    gmod.set_param_hint('x1', value=bgParsPos[1])
                    gmod.set_param_hint('y1', value=bgParsPos[2])
                    if separateNegParams:
                        gmod.set_param_hint('x1Neg', value=bgParsNeg[1])
                        gmod.set_param_hint('y1Neg', value=bgParsNeg[2])
                if bgGradientOrder >= 2:
                    gmod.set_param_hint('xy', value=bgParsPos[3])
                    gmod.set_param_hint('x2', value=bgParsPos[4])
                    gmod.set_param_hint('y2', value=bgParsPos[5])
                    if separateNegParams:
                        gmod.set_param_hint('xyNeg', value=bgParsNeg[3])
                        gmod.set_param_hint('x2Neg', value=bgParsNeg[4])
                        gmod.set_param_hint('y2Neg', value=bgParsNeg[5])

        y, x = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
        in_x = np.array([x, y]).astype(np.float64)
        in_x[0, :] -= in_x[0, :].mean()  # center it!
        in_x[1, :] -= in_x[1, :].mean()

        # Instead of explicitly using a mask to ignore flagged pixels, just set the ignored pixels'
        #  weights to 0 in the fit. TBD: need to inspect mask planes to set this mask.
        mask = np.ones_like(z, dtype=bool)  # TBD: set mask values to False if the pixels are to be ignored

        # I'm not sure about the variance planes in the diffim (or convolved pre-sub. images
        # for that matter) so for now, let's just do an un-weighted least-squares fit
        # (override weights computed above).
        weights = mask.astype(np.float64)
        if self.posImage is not None and rel_weight > 0.:
            weights = np.array([np.ones_like(diArr), np.ones_like(diArr)*rel_weight,
                                np.ones_like(diArr)*rel_weight])

        # Set the weights to zero if mask is False
        if np.any(~mask):
            weights[~mask] = 0.

        # Filter out any nans, and make the weights 0.
        nans = (np.isnan(z) | np.isnan(weights))
        nNans = nans.sum()
        if nNans > 0:
            if nNans < len(z):
                z[nans] = np.nanmedian(z)
            else:
                z[nans] = 0
        weights[nans] = 0

        # Note that although we can, we're not required to set initial values for params here,
        # since we set their param_hint's above.
        # Can add "method" param to not use 'leastsq' (==levenberg-marquardt), e.g. "method='nelder'"
        with warnings.catch_warnings():
            # Ignore lmfit unknown argument warnings:
            # "psf, rel_weight, footprint, modelObj" all become pass-through kwargs for makeModel.
            warnings.filterwarnings("ignore", "The keyword argument .* does not match", UserWarning)
            result = gmod.fit(z, weights=weights, x=in_x, max_nfev=250,
                              method="leastsq",  # TODO: try using `least_squares` here for speed/robustness
                              verbose=verbose,
                              # see scipy docs for the meaning of these keywords
                              fit_kws={'ftol': tol, 'xtol': tol, 'gtol': tol,
                                       # Our model is float32 internally, so we need a larger epsfcn.
                                       'epsfcn': 1e-8},
                              psf=self.diffim.getPsf(),  # hereon: kwargs that get passed to makeModel()
                              rel_weight=rel_weight,
                              footprint=fp,
                              modelObj=dipoleModel)

        if verbose:  # the ci_report() seems to fail if neg params are constrained -- TBD why.
            # Never wanted in production - this takes a long time (longer than the fit!)
            # This is how to get confidence intervals out:
            #    https://lmfit.github.io/lmfit-py/confidence.html and
            #    http://cars9.uchicago.edu/software/python/lmfit/model.html
            print(result.fit_report(show_correl=False))
            if separateNegParams:
                print(result.ci_report())

        return result

    def fitDipole(self, source, tol=1e-7, rel_weight=0.1,
                  fitBackground=1, maxSepInSigma=5., separateNegParams=True,
                  bgGradientOrder=1, verbose=False, display=False):
        """Fit a dipole model to an input ``diaSource`` (wraps `fitDipoleImpl`).

        Actually, fits the subimage bounded by the input source's
        footprint) and optionally constrain the fit using the
        pre-subtraction images self.posImage (science) and
        self.negImage (template). Wraps the output into a
        `pipeBase.Struct` named tuple after computing additional
        statistics such as orientation and SNR.

        Parameters
        ----------
        source : `lsst.afw.table.SourceRecord`
            Record containing the (merged) dipole source footprint detected on the diffim
        tol : `float`, optional
            Tolerance parameter for scipy.leastsq() optimization
        rel_weight : `float`, optional
            Weighting of posImage/negImage relative to the diffim in the fit
        fitBackground : `int`, {0, 1, 2}, optional
            How to fit linear background gradient in posImage/negImage

                - 0: do not fit background at all
                - 1 (default): pre-fit the background using linear least squares and then do not fit it
                  as part of the dipole fitting optimization
                - 2: pre-fit the background using linear least squares (as in 1), and use the parameter
                  estimates from that fit as starting parameters for an integrated "re-fit" of the
                  background as part of the overall dipole fitting optimization.
        maxSepInSigma : `float`, optional
            Allowed window of centroid parameters relative to peak in input source footprint
        separateNegParams : `bool`, optional
            Fit separate parameters to the flux and background gradient in
        bgGradientOrder : `int`, {0, 1, 2}, optional
            Desired polynomial order of background gradient
        verbose: `bool`, optional
            Be verbose
        display
            Display input data, best fit model(s) and residuals in a matplotlib window.

        Returns
        -------
        result : `struct`
            `pipeBase.Struct` object containing the fit parameters and other information.

        result : `callable`
            `lmfit.MinimizerResult` object for debugging and error estimation, etc.

        Notes
        -----
        Parameter `fitBackground` has three options, thus it is an integer:

        """

        fitResult = self.fitDipoleImpl(
            source, tol=tol, rel_weight=rel_weight, fitBackground=fitBackground,
            maxSepInSigma=maxSepInSigma, separateNegParams=separateNegParams,
            bgGradientOrder=bgGradientOrder, verbose=verbose)

        # Display images, model fits and residuals (currently uses matplotlib display functions)
        if display:
            fp = source.getFootprint()
            self.displayFitResults(fp, fitResult)

        # usually around 0.1 -- the minimum flux allowed -- i.e. bad fit.
        if fitResult.params['flux'].value <= 1.:
            return None, fitResult

        # TODO: We could include covariances, which could be derived from
        # `fitResult.params[name].correl`, but those are correlations.
        posCentroid = measBase.CentroidResult(fitResult.params['xcenPos'].value,
                                              fitResult.params['ycenPos'].value,
                                              fitResult.params['xcenPos'].stderr,
                                              fitResult.params['ycenPos'].stderr)
        negCentroid = measBase.CentroidResult(fitResult.params['xcenNeg'].value,
                                              fitResult.params['ycenNeg'].value,
                                              fitResult.params['xcenNeg'].stderr,
                                              fitResult.params['ycenNeg'].stderr)
        centroid = measBase.CentroidResult((fitResult.params['xcenPos'] + fitResult.params['xcenNeg']) / 2,
                                           (fitResult.params['ycenPos'] + fitResult.params['ycenNeg']) / 2.,
                                           math.sqrt(posCentroid.xErr**2 + negCentroid.xErr**2),
                                           math.sqrt(posCentroid.yErr**2 + negCentroid.yErr**2))
        dx = fitResult.params['xcenPos'].value - fitResult.params['xcenNeg'].value
        dy = fitResult.params['ycenPos'].value - fitResult.params['ycenNeg'].value
        angle = np.arctan2(dy, dx)

        # Exctract flux value, compute signalToNoise from flux/variance_within_footprint
        # Also extract the stderr of flux estimate.
        def computeSumVariance(exposure, footprint):
            return math.sqrt(np.nansum(exposure[footprint.getBBox(), afwImage.PARENT].variance.array))

        fluxVal = fluxVar = fitParams['flux']
        fluxVal = fluxVar = fitResult.params['flux'].value
        fluxErr = fluxErrNeg = fitResult.params['flux'].stderr
        if self.posImage is not None:
            fluxVar = computeSumVariance(self.posImage, source.getFootprint())
        else:
            fluxVar = computeSumVariance(self.diffim, source.getFootprint())

        fluxValNeg, fluxVarNeg = fluxVal, fluxVar
        if separateNegParams:
            fluxValNeg = fitResult.params['fluxNeg'].value
            fluxErrNeg = fitResult.params['fluxNeg'].stderr
        if self.negImage is not None:
            fluxVarNeg = computeSumVariance(self.negImage, source.getFootprint())

        try:
            signalToNoise = math.sqrt((fluxVal/fluxVar)**2 + (fluxValNeg/fluxVarNeg)**2)
        except ZeroDivisionError:  # catch divide by zero - should never happen.
            signalToNoise = np.nan

        out = Struct(posCentroid=posCentroid, negCentroid=negCentroid, centroid=centroid,
                     posFlux=fluxVal, negFlux=-fluxValNeg, posFluxErr=fluxErr, negFluxErr=fluxErrNeg,
                     orientation=angle,
                     signalToNoise=signalToNoise, chi2=fitResult.chisqr, redChi2=fitResult.redchi,
                     nData=fitResult.ndata)

        # fitResult may be returned for debugging
        return out, fitResult

    def displayFitResults(self, footprint, result):
        """Display data, model fits and residuals (currently uses matplotlib display functions).

        Parameters
        ----------
        footprint : `lsst.afw.detection.Footprint`
            Footprint containing the dipole that was fit
        result : `lmfit.MinimizerResult`
            `lmfit.MinimizerResult` object returned by `lmfit` optimizer
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as err:
            self.log.warning('Unable to import matplotlib: %s', err)
            raise err

        def display2dArray(ax, arr, x, y, xErr, yErr, title, extent=None):
            """Use `matplotlib.pyplot.imshow` to display a 2-D array with a given coordinate range.
            """
            fig = ax.imshow(arr, origin='lower', interpolation='none', cmap='gray', extent=extent)
            ax.set_title(title)
            ax.errorbar(x["total"], y["total"], xErr["total"], yErr["total"], c="cyan")
            ax.errorbar(x["Pos"], y["Pos"], xErr["Pos"], yErr["Pos"], c="green")
            ax.errorbar(x["Neg"], y["Neg"], xErr["Neg"], yErr["Neg"], c="red")
            return fig

        z = result.data
        fit = result.best_fit
        bbox = footprint.getBBox()
        extent = (bbox.getBeginX(), bbox.getEndX(), bbox.getBeginY(), bbox.getEndY())

        if z.shape[0] == 3:
            x, y, xErr, yErr = {}, {}, {}, {}
            for name in ("Pos", "Neg"):
                x[name] = result.best_values[f"xcen{name}"]
                y[name] = result.best_values[f"ycen{name}"]
                xErr[name] = result.params[f"xcen{name}"].stderr
                yErr[name] = result.params[f"ycen{name}"].stderr
            x["total"] = (result.best_values["xcenPos"] + result.best_values["xcenNeg"])/2
            y["total"] = (result.best_values["ycenPos"] + result.best_values["ycenNeg"])/2
            xErr["total"] = math.sqrt(xErr["Pos"]**2 + xErr["Neg"]**2)
            yErr["total"] = math.sqrt(yErr["Pos"]**2 + yErr["Neg"]**2)

            fig, axes = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(8, 8))
            for i, label in enumerate(("total", "Pos", "Neg")):
                display2dArray(axes[i][0], z[i, :], x, y, xErr, yErr,
                               f'Data {label}', extent=extent)
                display2dArray(axes[i][1], fit[i, :], x, y, xErr, yErr,
                               f'Model {label}', extent=extent)
                display2dArray(axes[i][2], z[i, :] - fit[i, :], x, y, xErr, yErr,
                               f'Residual {label}', extent=extent)

                plt.setp(axes[i][1].get_yticklabels(), visible=False)
                plt.setp(axes[i][2].get_yticklabels(), visible=False)
                if i != 2:  # remove top two row x-axis labels
                    plt.setp(axes[i][0].get_xticklabels(), visible=False)
                    plt.setp(axes[i][1].get_xticklabels(), visible=False)
                    plt.setp(axes[i][2].get_xticklabels(), visible=False)
        else:
            fig, axes = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(10, 2.5))
            display2dArray(axes[0], z, 'Data', extent=extent)
            display2dArray(axes[1], 'Model', extent=extent)
            display2dArray(axes[2], z - fit, 'Residual', extent=extent)

        fig.tight_layout(pad=0, w_pad=0, h_pad=0)
        plt.show()


@measBase.register("ip_diffim_DipoleFit")
class DipoleFitPlugin(measBase.SingleFramePlugin):
    """A single frame measurement plugin that fits dipoles to all merged (two-peak) ``diaSources``.

    This measurement plugin accepts up to three input images in
    its `measure` method. If these are provided, it includes data
    from the pre-subtraction posImage (science image) and optionally
    negImage (template image) to constrain the fit. The meat of the
    fitting routines are in the class `~lsst.module.name.DipoleFitAlgorithm`.

    Notes
    -----
    The motivation behind this plugin and the necessity for including more than
    one exposure are documented in DMTN-007 (http://dmtn-007.lsst.io).

    This class is named `ip_diffim_DipoleFit` so that it may be used alongside
    the existing `ip_diffim_DipoleMeasurement` classes until such a time as those
    are deemed to be replaceable by this.
    """

    ConfigClass = DipoleFitPluginConfig
    DipoleFitAlgorithmClass = DipoleFitAlgorithm  # Pointer to the class that performs the fit

    FAILURE_EDGE = 1   # too close to the edge
    FAILURE_FIT = 2    # failure in the fitting
    FAILURE_NOT_DIPOLE = 4  # input source is not a putative dipole to begin with
    FAILURE_TOO_LARGE = 8  # input source is too large to be fit

    @classmethod
    def getExecutionOrder(cls):
        """This algorithm simultaneously fits the centroid and flux, and does
        not require any previous centroid fit.
        """
        return cls.CENTROID_ORDER

    def __init__(self, config, name, schema, metadata, logName=None):
        if logName is None:
            logName = name
        measBase.SingleFramePlugin.__init__(self, config, name, schema, metadata, logName=logName)

        self.log = logging.getLogger(logName)

        self._setupSchema(config, name, schema, metadata)

    def _setupSchema(self, config, name, schema, metadata):
        """Add fields for the outputs, and save the keys for fast assignment.
        """
        self.posFluxKey = measBase.FluxResultKey.addFields(schema,
                                                           schema.join(name, "pos"),
                                                           "Dipole positive lobe instrumental flux.")
        self.negFluxKey = measBase.FluxResultKey.addFields(schema,
                                                           schema.join(name, "neg"),
                                                           "Dipole negative lobe instrumental flux.")
        doc = "Dipole overall instrumental flux (mean of absolute value of positive and negative lobes)."
        self.fluxKey = measBase.FluxResultKey.addFields(schema, name, doc)

        self.posCentroidKey = measBase.CentroidResultKey.addFields(schema,
                                                                   schema.join(name, "pos"),
                                                                   "Dipole positive lobe centroid position.",
                                                                   measBase.UncertaintyEnum.SIGMA_ONLY)
        self.negCentroidKey = measBase.CentroidResultKey.addFields(schema,
                                                                   schema.join(name, "neg"),
                                                                   "Dipole negative lobe centroid position.",
                                                                   measBase.UncertaintyEnum.SIGMA_ONLY)
        self.centroidKey = measBase.CentroidResultKey.addFields(schema,
                                                                name,
                                                                "Dipole centroid position.",
                                                                measBase.UncertaintyEnum.SIGMA_ONLY)

        self.orientationKey = schema.addField(
            schema.join(name, "orientation"), type=float, units="rad",
            doc="Dipole orientation.  Convention is CCW from +x on image.")

        self.separationKey = schema.addField(
            schema.join(name, "separation"), type=float, units="pixel",
            doc="Pixel separation between positive and negative lobes of dipole")

        self.chi2dofKey = schema.addField(
            schema.join(name, "chi2dof"), type=float,
            doc="Chi2 per degree of freedom (chi2/(nData-nVariables)) of dipole fit")

        self.nDataKey = schema.addField(
            schema.join(name, "nData"), type=np.int64,
            doc="Number of data points in the dipole fit")

        self.signalToNoiseKey = schema.addField(
            schema.join(name, "signalToNoise"), type=float,
            doc="Estimated signal-to-noise of dipole fit")

        self.classificationFlagKey = schema.addField(
            schema.join(name, "classification"), type="Flag",
            doc="Flag indicating diaSource is classified as a dipole")

        self.classificationAttemptedFlagKey = schema.addField(
            schema.join(name, "classificationAttempted"), type="Flag",
            doc="Flag indicating diaSource was attempted to be classified as a dipole")

        self.flagKey = schema.addField(
            schema.join(name, "flag"), type="Flag",
            doc="General failure flag for dipole fit")

        self.edgeFlagKey = schema.addField(
            schema.join(name, "flag", "edge"), type="Flag",
            doc="Flag set when dipole is too close to edge of image")

    def measureDipoles(self, measRecord, exposure, posExp=None, negExp=None):
        """Perform the non-linear least squares minimization on the putative dipole source.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            diaSources that will be measured using dipole measurement
        exposure : `lsst.afw.image.Exposure`
            Difference exposure on which the diaSources were detected; `exposure = posExp-negExp`
            If both `posExp` and `negExp` are `None`, will attempt to fit the
            dipole to just the `exposure` with no constraint.
        posExp : `lsst.afw.image.Exposure`, optional
            "Positive" exposure, typically a science exposure, or None if unavailable
            When `posExp` is `None`, will compute `posImage = exposure + negExp`.
        negExp : `lsst.afw.image.Exposure`, optional
            "Negative" exposure, typically a template exposure, or None if unavailable
            When `negExp` is `None`, will compute `negImage = posExp - exposure`.

        Notes
        -----
        The main functionality of this routine was placed outside of
        this plugin (into `DipoleFitAlgorithm.fitDipole()`) so that
        `DipoleFitAlgorithm.fitDipole()` can be called separately for
        testing (@see `tests/testDipoleFitter.py`)
        """
        result = None
        pks = measRecord.getFootprint().getPeaks()

        # Check if the footprint consists of a putative dipole - else don't fit it.
        if (
                # One peak in the footprint (not a dipole)
                ((nPeaks := len(pks)) <= 1)
                # Peaks are the same sign (not a dipole); peaks are ordered
                # from highest to lowest.
                or (nPeaks > 1 and (np.sign(pks[0].getPeakValue())
                    == np.sign(pks[-1].getPeakValue())))
        ):
            if not self.config.fitAllDiaSources:
                # Non-dipoles fall back on the centroid slot for positions,
                # errors, and the failure flag, if we're not fitting them.
                measRecord[self.centroidKey.getX()] = measRecord.getX()
                measRecord[self.centroidKey.getY()] = measRecord.getY()
                self.centroidKey.getCentroidErr().set(measRecord, measRecord.getCentroidErr())
                measRecord[self.flagKey] = measRecord.getCentroidFlag()
                return

        # Footprint is too large (not a dipole).
        if ((area := measRecord.getFootprint().getArea()) > self.config.maxFootprintArea):
            self.fail(measRecord, measBase.MeasurementError(f"{area} > {self.config.maxFootprintArea}",
                                                            self.FAILURE_TOO_LARGE))
            return

        try:
            alg = self.DipoleFitAlgorithmClass(exposure, posImage=posExp, negImage=negExp)
            result, _ = alg.fitDipole(
                measRecord, rel_weight=self.config.relWeight,
                tol=self.config.tolerance,
                maxSepInSigma=self.config.maxSeparation,
                fitBackground=self.config.fitBackground,
                separateNegParams=self.config.fitSeparateNegParams,
                verbose=False, display=False)
        except pexExcept.LengthError:
            self.fail(measRecord, measBase.MeasurementError('edge failure', self.FAILURE_EDGE))
        except Exception as e:
            errorMessage = f"Exception in dipole fit. {e.__class__.__name__}: {e}"
            self.fail(measRecord, measBase.MeasurementError(errorMessage, self.FAILURE_FIT))

        self.log.debug("Dipole fit result: %d %s", measRecord.getId(), str(result))

        if result is None:
            self.fail(measRecord, measBase.MeasurementError("bad dipole fit", self.FAILURE_FIT))
            return

        # add chi2, coord/flux uncertainties (TBD), dipole classification
        # Add the relevant values to the measRecord
        measRecord[self.posFluxKey.getInstFlux()] = result.posFlux
        measRecord[self.posFluxKey.getInstFluxErr()] = result.signalToNoise   # to be changed to actual sigma!
        self.posCentroidKey.set(measRecord, result.posCentroid)
        measRecord[self.posCentroidKey.getY()] = result.posCentroidY

        measRecord[self.negFluxKey.getInstFlux()] = result.negFlux
        measRecord[self.negFluxKey.getInstFluxErr()] = result.signalToNoise   # to be changed to actual sigma!
        self.negCentroidKey.set(measRecord, result.negCentroid)
        measRecord[self.negCentroidKey.getY()] = result.negCentroidY

        # Dia source flux: average of pos+neg
        measRecord[self.fluxKey.getInstFlux()] = (abs(result.posFlux) + abs(result.negFlux))/2.
        measRecord[self.orientationKey] = result.orientation
        measRecord[self.separationKey] = math.sqrt((result.posCentroid.x - result.negCentroid.x)**2
                                                   + (result.posCentroid.y - result.negCentroid.y)**2)

        self.centroidKey.set(measRecord, result.centroid)

        measRecord[self.signalToNoiseKey] = result.signalToNoise
        measRecord[self.chi2dofKey] = result.redChi2

        if result.nData >= 1:
            measRecord[self.nDataKey] = result.nData
        else:
            measRecord[self.nDataKey] = 0

        self.doClassify(measRecord, result.chi2)

    def doClassify(self, measRecord, chi2val):
        """Classify a source as a dipole.

        Parameters
        ----------
        measRecord : TODO: DM-17458
            TODO: DM-17458
        chi2val : TODO: DM-17458
            TODO: DM-17458

        Notes
        -----
        Sources are classified as dipoles, or not, according to three criteria:

        1. Does the total signal-to-noise surpass the ``minSn``?
        2. Are the pos/neg fluxes greater than 1.0 and no more than 0.65 (``maxFluxRatio``)
           of the total flux? By default this will never happen since ``posFlux == negFlux``.
        3. Is it a good fit (``chi2dof`` < 1)? (Currently not used.)
        """

        # First, does the total signal-to-noise surpass the minSn?
        passesSn = measRecord[self.signalToNoiseKey] > self.config.minSn

        # Second, are the pos/neg fluxes greater than 1.0 and no more than 0.65 (param maxFluxRatio)
        # of the total flux? By default this will never happen since posFlux = negFlux.
        passesFluxPos = (abs(measRecord[self.posFluxKey.getInstFlux()])
                         / (measRecord[self.fluxKey.getInstFlux()]*2.)) < self.config.maxFluxRatio
        passesFluxPos &= (abs(measRecord[self.posFluxKey.getInstFlux()]) >= 1.0)
        passesFluxNeg = (abs(measRecord[self.negFluxKey.getInstFlux()])
                         / (measRecord[self.fluxKey.getInstFlux()]*2.)) < self.config.maxFluxRatio
        passesFluxNeg &= (abs(measRecord[self.negFluxKey.getInstFlux()]) >= 1.0)
        allPass = (passesSn and passesFluxPos and passesFluxNeg)  # and passesChi2)

        # Third, is it a good fit (chi2dof < 1)?
        # Use scipy's chi2 cumulative distrib to estimate significance
        # This doesn't really work since I don't trust the values in the variance plane (which
        #   affects the least-sq weights, which affects the resulting chi2).
        # But I'm going to keep this here for future use.
        if False:
            from scipy.stats import chi2
            ndof = chi2val / measRecord[self.chi2dofKey]
            significance = chi2.cdf(chi2val, ndof)
            passesChi2 = significance < self.config.maxChi2DoF
            allPass = allPass and passesChi2

        measRecord.set(self.classificationAttemptedFlagKey, True)

        if allPass:  # Note cannot pass `allPass` into the `measRecord.set()` call below...?
            measRecord.set(self.classificationFlagKey, True)
        else:
            measRecord.set(self.classificationFlagKey, False)

    def fail(self, measRecord, error=None):
        """Catch failures and set the correct flags.

        Fallback on the current slot centroid positions, but set the dipole
        failure flag, since we attempted to fit the source.
        """
        measRecord[self.centroidKey.getX()] = measRecord.getX()
        measRecord[self.centroidKey.getY()] = measRecord.getY()
        self.centroidKey.getCentroidErr().set(measRecord, measRecord.getCentroidErr())

        measRecord.set(self.flagKey, True)
        if error is not None:
            if error.getFlagBit() == self.FAILURE_EDGE:
                self.log.debug('DipoleFitPlugin not run on record %d: %s', measRecord.getId(), str(error))
                measRecord.set(self.edgeFlagKey, True)
            if error.getFlagBit() == self.FAILURE_FIT:
                self.log.warning('DipoleFitPlugin failed on record %d: %s', measRecord.getId(), str(error))
            if error.getFlagBit() == self.FAILURE_TOO_LARGE:
                self.log.debug('DipoleFitPlugin not run on record with too large footprint %d: %s',
                               measRecord.getId(), str(error))
        else:
            self.log.warning('DipoleFitPlugin failed on record %d', measRecord.getId())
