# This file is part of ip_diffim.
#
# LSST Data Management System
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
# See COPYRIGHT file at the top of the source tree.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

import numpy as np
from scipy import ndimage
from lsst.afw.coord import differentialRefraction
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.pipe.base as pipeBase

__all__ = ["DcrModel", "applyDcr", "calculateDcr", "calculateImageParallacticAngle"]


class DcrModel:
    """A model of the true sky after correcting chromatic effects.

    Attributes
    ----------
    dcrNumSubfilters : `int`
        Number of sub-filters used to model chromatic effects within a band.
    modelImages : `list` of `lsst.afw.image.Image`
        A list of masked images, each containing the model for one subfilter

    Notes
    -----
    The ``DcrModel`` contains an estimate of the true sky, at a higher
    wavelength resolution than the input observations. It can be forward-
    modeled to produce Differential Chromatic Refraction (DCR) matched
    templates for a given ``Exposure``, and provides utilities for conditioning
    the model in ``dcrAssembleCoadd`` to avoid oscillating solutions between
    iterations of forward modeling or between the subfilters of the model.
    """

    def __init__(self, modelImages, effectiveWavelength, bandwidth, filterLabel=None, psf=None,
                 bbox=None, wcs=None, mask=None, variance=None, photoCalib=None):
        self.dcrNumSubfilters = len(modelImages)
        self.modelImages = modelImages
        self._filterLabel = filterLabel
        self._effectiveWavelength = effectiveWavelength
        self._bandwidth = bandwidth
        self._psf = psf
        self._bbox = bbox
        self._wcs = wcs
        self._mask = mask
        self._variance = variance
        self.photoCalib = photoCalib

    @classmethod
    def fromImage(cls, maskedImage, dcrNumSubfilters, effectiveWavelength, bandwidth,
                  wcs=None, filterLabel=None, psf=None, photoCalib=None):
        """Initialize a DcrModel by dividing a coadd between the subfilters.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            Input coadded image to divide equally between the subfilters.
        dcrNumSubfilters : `int`
            Number of sub-filters used to model chromatic effects within a
            band.
        effectiveWavelength : `float`
            The effective wavelengths of the current filter, in nanometers.
        bandwidth : `float`
            The bandwidth of the current filter, in nanometers.
        wcs : `lsst.afw.geom.SkyWcs`
            Coordinate system definition (wcs) for the exposure.
        filterLabel : `lsst.afw.image.FilterLabel`, optional
            The filter label, set in the current instruments' obs package.
            Required for any calculation of DCR, including making matched
            templates.
        psf : `lsst.afw.detection.Psf`, optional
            Point spread function (PSF) of the model.
            Required if the ``DcrModel`` will be persisted.
        photoCalib : `lsst.afw.image.PhotoCalib`, optional
            Calibration to convert instrumental flux and
            flux error to nanoJansky.

        Returns
        -------
        dcrModel : `lsst.pipe.tasks.DcrModel`
            Best fit model of the true sky after correcting chromatic effects.
        """
        # NANs will potentially contaminate the entire image,
        # depending on the shift or convolution type used.
        model = maskedImage.image.clone()
        mask = maskedImage.mask.clone()
        bbox = maskedImage.getBBox()
        # We divide the variance by N and not N**2 because we will assume each
        # subfilter is independent. That means that the significance of
        # detected sources will be lower by a factor of sqrt(N) in the
        # subfilter images, but we will recover it when we combine the
        # subfilter images to construct matched templates.
        variance = maskedImage.variance.clone()
        variance /= dcrNumSubfilters
        model /= dcrNumSubfilters
        modelImages = [model, ]
        for subfilter in range(1, dcrNumSubfilters):
            modelImages.append(model.clone())
        return cls(modelImages, effectiveWavelength, bandwidth,
                   filterLabel=filterLabel, psf=psf, bbox=bbox, wcs=wcs,
                   mask=mask, variance=variance, photoCalib=photoCalib)

    @classmethod
    def fromQuantum(cls, availableCoaddRefs, effectiveWavelength, bandwidth, numSubfilters):
        """Load an existing DcrModel from a Gen 3 repository.

        Parameters
        ----------
        availableCoaddRefs : `dict` [`int`, `lsst.daf.butler.DeferredDatasetHandle`]
            Dictionary of spatially relevant retrieved coadd patches,
            indexed by their sequential patch number.
        effectiveWavelength : `float`
            The effective wavelengths of the current filter, in nanometers.
        bandwidth : `float`
            The bandwidth of the current filter, in nanometers.
        numSubfilters : `int`
            Number of subfilters in the DcrCoadd.

        Returns
        -------
        dcrModel : `lsst.pipe.tasks.DcrModel`
            Best fit model of the true sky after correcting chromatic effects.
        """
        filterLabel = None
        psf = None
        bbox = None
        wcs = None
        mask = None
        variance = None
        photoCalib = None
        modelImages = [None]*numSubfilters

        for coaddRef in availableCoaddRefs:
            subfilter = coaddRef.dataId["subfilter"]
            dcrCoadd = coaddRef.get()
            if filterLabel is None:
                filterLabel = dcrCoadd.getFilter()
            if psf is None:
                psf = dcrCoadd.getPsf()
            if bbox is None:
                bbox = dcrCoadd.getBBox()
            if wcs is None:
                wcs = dcrCoadd.wcs
            if mask is None:
                mask = dcrCoadd.mask
            if variance is None:
                variance = dcrCoadd.variance
            if photoCalib is None:
                photoCalib = dcrCoadd.getPhotoCalib()
            modelImages[subfilter] = dcrCoadd.image
        return cls(modelImages, effectiveWavelength, bandwidth, filterLabel,
                   psf, bbox, wcs, mask, variance, photoCalib)

    def __len__(self):
        """Return the number of subfilters.

        Returns
        -------
        dcrNumSubfilters : `int`
            The number of DCR subfilters in the model.
        """
        return self.dcrNumSubfilters

    def __getitem__(self, subfilter):
        """Iterate over the subfilters of the DCR model.

        Parameters
        ----------
        subfilter : `int`
            Index of the current ``subfilter`` within the full band.
            Negative indices are allowed, and count in reverse order
            from the highest ``subfilter``.

        Returns
        -------
        modelImage : `lsst.afw.image.Image`
            The DCR model for the given ``subfilter``.

        Raises
        ------
        IndexError
            If the requested ``subfilter`` is greater or equal to the number
            of subfilters in the model.
        """
        if np.abs(subfilter) >= len(self):
            raise IndexError("subfilter out of bounds.")
        return self.modelImages[subfilter]

    def __setitem__(self, subfilter, maskedImage):
        """Update the model image for one subfilter.

        Parameters
        ----------
        subfilter : `int`
            Index of the current subfilter within the full band.
        maskedImage : `lsst.afw.image.Image`
            The DCR model to set for the given ``subfilter``.

        Raises
        ------
        IndexError
            If the requested ``subfilter`` is greater or equal to the number
            of subfilters in the model.
        ValueError
            If the bounding box of the new image does not match.
        """
        if np.abs(subfilter) >= len(self):
            raise IndexError("subfilter out of bounds.")
        if maskedImage.getBBox() != self.bbox:
            raise ValueError("The bounding box of a subfilter must not change.")
        self.modelImages[subfilter] = maskedImage

    @property
    def effectiveWavelength(self):
        """Return the effective wavelength of the model.

        Returns
        -------
        effectiveWavelength : `float`
            The effective wavelength of the current filter, in nanometers.
        """
        return self._effectiveWavelength

    @property
    def filter(self):
        """Return the filter label for the model.

        Returns
        -------
        filterLabel : `lsst.afw.image.FilterLabel`
            The filter used for the input observations.
        """
        return self._filterLabel

    @property
    def bandwidth(self):
        """Return the bandwidth of the model.

        Returns
        -------
        bandwidth : `float`
            The bandwidth of the current filter, in nanometers.
        """
        return self._bandwidth

    @property
    def psf(self):
        """Return the psf of the model.

        Returns
        -------
        psf : `lsst.afw.detection.Psf`
            Point spread function (PSF) of the model.
        """
        return self._psf

    @property
    def bbox(self):
        """Return the common bounding box of each subfilter image.

        Returns
        -------
        bbox : `lsst.afw.geom.Box2I`
            Bounding box of the DCR model.
        """
        return self._bbox

    @property
    def wcs(self):
        """Return the WCS of each subfilter image.

        Returns
        -------
        wcs : `lsst.afw.geom.SkyWcs`
            Coordinate system definition (wcs) for the exposure.
        """
        return self._wcs

    @property
    def mask(self):
        """Return the common mask of each subfilter image.

        Returns
        -------
        mask : `lsst.afw.image.Mask`
            Mask plane of the DCR model.
        """
        return self._mask

    @property
    def variance(self):
        """Return the common variance of each subfilter image.

        Returns
        -------
        variance : `lsst.afw.image.Image`
            Variance plane of the DCR model.
        """
        return self._variance

    def getReferenceImage(self, bbox=None):
        """Calculate a reference image from the average of the subfilter
        images.

        Parameters
        ----------
        bbox : `lsst.afw.geom.Box2I`, optional
            Sub-region of the coadd. Returns the entire image if `None`.

        Returns
        -------
        refImage : `numpy.ndarray`
            The reference image with no chromatic effects applied.
        """
        bbox = bbox or self.bbox
        return np.mean([model[bbox].array for model in self], axis=0)

    def assign(self, dcrSubModel, bbox=None):
        """Update a sub-region of the ``DcrModel`` with new values.

        Parameters
        ----------
        dcrSubModel : `lsst.pipe.tasks.DcrModel`
            New model of the true scene after correcting chromatic effects.
        bbox : `lsst.afw.geom.Box2I`, optional
            Sub-region of the coadd.
            Defaults to the bounding box of ``dcrSubModel``.

        Raises
        ------
        ValueError
            If the new model has a different number of subfilters.
        """
        if len(dcrSubModel) != len(self):
            raise ValueError("The number of DCR subfilters must be the same "
                             "between the old and new models.")
        bbox = bbox or self.bbox
        for model, subModel in zip(self, dcrSubModel):
            model.assign(subModel[bbox], bbox)

    def buildMatchedTemplate(self, exposure=None, order=3,
                             visitInfo=None, bbox=None, mask=None,
                             splitSubfilters=True, splitThreshold=0., amplifyModel=1.):
        """Create a DCR-matched template image for an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`, optional
            The input exposure to build a matched template for.
            May be omitted if all of the metadata is supplied separately
        order : `int`, optional
            Interpolation order of the DCR shift.
        visitInfo : `lsst.afw.image.VisitInfo`, optional
            Metadata for the exposure. Ignored if ``exposure`` is set.
        bbox : `lsst.afw.geom.Box2I`, optional
            Sub-region of the coadd, or use the entire coadd if not supplied.
        mask : `lsst.afw.image.Mask`, optional
            reference mask to use for the template image.
        splitSubfilters : `bool`, optional
            Calculate DCR for two evenly-spaced wavelengths in each subfilter,
            instead of at the midpoint. Default: True
        splitThreshold : `float`, optional
            Minimum DCR difference within a subfilter required to use
            ``splitSubfilters``
        amplifyModel : `float`, optional
            Multiplication factor to amplify differences between model planes.
            Used to speed convergence of iterative forward modeling.

        Returns
        -------
        templateImage : `lsst.afw.image.ImageF`
            The DCR-matched template

        Raises
        ------
        ValueError
            If neither ``exposure`` or ``visitInfo`` are set.
        """
        if self.effectiveWavelength is None or self.bandwidth is None:
            raise ValueError("'effectiveWavelength' and 'bandwidth' must be set for the DcrModel in order "
                             "to calculate DCR.")
        if exposure is not None:
            visitInfo = exposure.getInfo().getVisitInfo()
        elif visitInfo is None:
            raise ValueError("Either exposure or visitInfo must be set.")
        if bbox is None:
            bbox = self.bbox
        else:
            bbox = bbox.clippedTo(self.bbox)
        dcrShift = calculateDcr(visitInfo, self.wcs, self.effectiveWavelength, self.bandwidth, len(self),
                                splitSubfilters=splitSubfilters, bbox=bbox)
        templateImage = afwImage.ImageF(bbox)
        refModel = None
        for subfilter, dcr in enumerate(dcrShift):
            if self[subfilter] is None:
                # It is possible to load only a single DCR subfilter at a time.
                self.log.debug("Skipping missing DCR model subfilter %d", subfilter)
                continue
            if amplifyModel > 1:
                if refModel is None:
                    # amplifyModel is only an option while constructing the DcrModel,
                    # and we don't want to calculate a reference image during image differencing.
                    refModel = self.getReferenceImage(bbox)
                model = (self[subfilter][bbox].array - refModel)*amplifyModel + refModel
            else:
                model = self[subfilter][bbox].array
            templateImage.array += applyDcr(model, dcr, splitSubfilters=splitSubfilters,
                                            splitThreshold=splitThreshold, order=order)
        return templateImage

    def buildMatchedExposure(self, exposure=None,
                             visitInfo=None, bbox=None, mask=None):
        """Wrapper to create an exposure from a template image.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`, optional
            The input exposure to build a matched template for.
            May be omitted if all of the metadata is supplied separately
        visitInfo : `lsst.afw.image.VisitInfo`, optional
            Metadata for the exposure. Ignored if ``exposure`` is set.
        bbox : `lsst.afw.geom.Box2I`, optional
            Sub-region of the coadd, or use the entire coadd if not supplied.
        mask : `lsst.afw.image.Mask`, optional
            reference mask to use for the template image.

        Returns
        -------
        templateExposure : `lsst.afw.image.exposureF`
            The DCR-matched template

        Raises
        ------
        RuntimeError
            If no `photcCalib` is set.
        """
        if bbox is None:
            bbox = self.bbox
        templateImage = self.buildMatchedTemplate(exposure=exposure, visitInfo=visitInfo,
                                                  bbox=bbox, mask=mask)
        maskedImage = afwImage.MaskedImageF(bbox)
        maskedImage.image = templateImage[bbox]
        maskedImage.mask = self.mask[bbox]
        maskedImage.variance = self.variance[bbox]
        # The variance of the stacked image will be `dcrNumSubfilters`
        # times the variance of the individual subfilters.
        maskedImage.variance *= self.dcrNumSubfilters
        templateExposure = afwImage.ExposureF(bbox, self.wcs)
        templateExposure.setMaskedImage(maskedImage[bbox])
        templateExposure.setPsf(self.psf)
        templateExposure.setFilter(self.filter)
        if self.photoCalib is None:
            raise RuntimeError("No PhotoCalib set for the DcrModel. "
                               "If the DcrModel was created from a masked image"
                               " you must also specify the photoCalib.")
        templateExposure.setPhotoCalib(self.photoCalib)
        return templateExposure

    def buildMatchedExposureHandle(self, exposure=None, visitInfo=None, bbox=None, mask=None):
        """Create in-memory butler dataset reference containing the DCR-matched
        template.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`, optional
            The input exposure to build a matched template for.
            May be omitted if all of the metadata is supplied separately
        visitInfo : `lsst.afw.image.VisitInfo`, optional
            Metadata for the exposure. Ignored if ``exposure`` is set.
        bbox : `lsst.afw.geom.Box2I`, optional
            Sub-region of the coadd, or use the entire coadd if not supplied.
        mask : `lsst.afw.image.Mask`, optional
            reference mask to use for the template image.

        Returns
        -------
        templateExposureHandle: `lsst.pipe.base.InMemoryDatasetHandle`
            In-memory butler dataset reference containing the DCR-matched
            template.

        Raises
        ------
        ValueError
            If no `exposure` or `visitInfo` is set.
        """
        if exposure is not None:
            visitInfo = exposure.visitInfo
        elif visitInfo is None:
            raise ValueError("Either exposure or visitInfo must be set.")
        templateExposure = self.buildMatchedExposure(
            exposure=exposure, visitInfo=visitInfo, bbox=bbox, mask=mask
        )
        templateExposureHandle = pipeBase.InMemoryDatasetHandle(
            templateExposure,
            storageClass="ExposureF",
            copy=False,
            photoCalib=self.photoCalib
        )
        return templateExposureHandle

    def conditionDcrModel(self, modelImages, bbox, gain=1.):
        """Average two iterations' solutions to reduce oscillations.

        Parameters
        ----------
        modelImages : `list` of `lsst.afw.image.Image`
            The new DCR model images from the current iteration.
            The values will be modified in place.
        bbox : `lsst.afw.geom.Box2I`
            Sub-region of the coadd
        gain : `float`, optional
            Relative weight to give the new solution when updating the model.
            Defaults to 1.0, which gives equal weight to both solutions.
        """
        # Calculate weighted averages of the images.
        for model, newModel in zip(self, modelImages):
            newModel *= gain
            newModel += model[bbox]
            newModel /= 1. + gain

    def regularizeModelIter(self, subfilter, newModel, bbox, regularizationFactor,
                            regularizationWidth=2):
        """Restrict large variations in the model between iterations.

        Parameters
        ----------
        subfilter : `int`
            Index of the current subfilter within the full band.
        newModel : `lsst.afw.image.Image`
            The new DCR model for one subfilter from the current iteration.
            Values in ``newModel`` that are extreme compared with the last
            iteration are modified in place.
        bbox : `lsst.afw.geom.Box2I`
            Sub-region to coadd
        regularizationFactor : `float`
            Maximum relative change of the model allowed between iterations.
        regularizationWidth : int, optional
            Minimum radius of a region to include in regularization, in pixels.
        """
        refImage = self[subfilter][bbox].array
        highThreshold = np.abs(refImage)*regularizationFactor
        lowThreshold = refImage/regularizationFactor
        newImage = newModel.array
        self.applyImageThresholds(newImage, highThreshold=highThreshold, lowThreshold=lowThreshold,
                                  regularizationWidth=regularizationWidth)

    def regularizeModelFreq(self, modelImages, bbox, statsCtrl, regularizationFactor,
                            regularizationWidth=2, mask=None, convergenceMaskPlanes="DETECTED"):
        """Restrict large variations in the model between subfilters.

        Parameters
        ----------
        modelImages : `list` of `lsst.afw.image.Image`
            The new DCR model images from the current iteration.
            The values will be modified in place.
        bbox : `lsst.afw.geom.Box2I`
            Sub-region to coadd
        statsCtrl : `lsst.afw.math.StatisticsControl`
            Statistics control object for coaddition.
        regularizationFactor : `float`
            Maximum relative change of the model allowed between subfilters.
        regularizationWidth : `int`, optional
            Minimum radius of a region to include in regularization, in pixels.
        mask : `lsst.afw.image.Mask`, optional
            Optional alternate mask
        convergenceMaskPlanes : `list` of `str`, or `str`, optional
            Mask planes to use to calculate convergence.

        Notes
        -----
        This implementation of frequency regularization restricts each
        subfilter image to be a smoothly-varying function times a reference
        image.
        """
        # ``regularizationFactor`` is the maximum change between subfilter
        # images, so the maximum difference between one subfilter image and the
        # average will be the square root of that.
        maxDiff = np.sqrt(regularizationFactor)
        noiseLevel = self.calculateNoiseCutoff(modelImages[0], statsCtrl, bufferSize=5, mask=mask, bbox=bbox)
        referenceImage = self.getReferenceImage(bbox)
        badPixels = np.isnan(referenceImage) | (referenceImage <= 0.)
        if np.sum(~badPixels) == 0:
            # Skip regularization if there are no valid pixels
            return
        referenceImage[badPixels] = 0.
        filterWidth = regularizationWidth
        fwhm = 2.*filterWidth
        # The noise should be lower in the smoothed image by
        # sqrt(Nsmooth) ~ fwhm pixels
        noiseLevel /= fwhm
        smoothRef = ndimage.gaussian_filter(referenceImage, filterWidth, mode='constant')
        # Add a three sigma offset to both the reference and model to prevent
        # dividing by zero. Note that this will also slightly suppress faint
        # variations in color.
        smoothRef += 3.*noiseLevel

        lowThreshold = smoothRef/maxDiff
        highThreshold = smoothRef*maxDiff
        for model in modelImages:
            self.applyImageThresholds(model.array,
                                      highThreshold=highThreshold,
                                      lowThreshold=lowThreshold,
                                      regularizationWidth=regularizationWidth)
            smoothModel = ndimage.gaussian_filter(model.array, filterWidth, mode='constant')
            smoothModel += 3.*noiseLevel
            relativeModel = smoothModel/smoothRef
            # Now sharpen the smoothed relativeModel using an alpha of 3.
            alpha = 3.
            relativeModel2 = ndimage.gaussian_filter(relativeModel, filterWidth/alpha)
            relativeModel += alpha*(relativeModel - relativeModel2)
            model.array = relativeModel*referenceImage

    def calculateNoiseCutoff(self, image, statsCtrl, bufferSize,
                             convergenceMaskPlanes="DETECTED", mask=None, bbox=None):
        """Helper function to calculate the background noise level of an image.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            The input image to evaluate the background noise properties.
        statsCtrl : `lsst.afw.math.StatisticsControl`
            Statistics control object for coaddition.
        bufferSize : `int`
            Number of additional pixels to exclude
            from the edges of the bounding box.
        convergenceMaskPlanes : `list` of `str`, or `str`
            Mask planes to use to calculate convergence.
        mask : `lsst.afw.image.Mask`, Optional
            Optional alternate mask
        bbox : `lsst.afw.geom.Box2I`, optional
            Sub-region of the masked image to calculate the noise level over.

        Returns
        -------
        noiseCutoff : `float`
            The threshold value to treat pixels as noise in an image..
        """
        if bbox is None:
            bbox = self.bbox
        if mask is None:
            mask = self.mask[bbox]
        bboxShrink = geom.Box2I(bbox)
        bboxShrink.grow(-bufferSize)
        convergeMask = mask.getPlaneBitMask(convergenceMaskPlanes)

        backgroundPixels = mask[bboxShrink].array & (statsCtrl.getAndMask() | convergeMask) == 0
        noiseCutoff = np.std(image[bboxShrink].array[backgroundPixels])
        return noiseCutoff

    def applyImageThresholds(self, image, highThreshold=None, lowThreshold=None, regularizationWidth=2):
        """Restrict image values to be between upper and lower limits.

        This method flags all pixels in an image that are outside of the given
        threshold values. The threshold values are taken from a reference
        image, so noisy pixels are likely to get flagged. In order to exclude
        those noisy pixels, the array of flags is eroded and dilated, which
        removes isolated pixels outside of the thresholds from the list of
        pixels to be modified. Pixels that remain flagged after this operation
        have their values set to the appropriate upper or lower threshold
        value.

        Parameters
        ----------
        image : `numpy.ndarray`
            The image to apply the thresholds to.
            The values will be modified in place.
        highThreshold : `numpy.ndarray`, optional
            Array of upper limit values for each pixel of ``image``.
        lowThreshold : `numpy.ndarray`, optional
            Array of lower limit values for each pixel of ``image``.
        regularizationWidth : `int`, optional
            Minimum radius of a region to include in regularization, in pixels.
        """
        # Generate the structure for binary erosion and dilation, which is used
        # to remove noise-like pixels. Groups of pixels with a radius smaller
        # than ``regularizationWidth`` will be excluded from regularization.
        filterStructure = ndimage.iterate_structure(ndimage.generate_binary_structure(2, 1),
                                                    regularizationWidth)
        if highThreshold is not None:
            highPixels = image > highThreshold
            if regularizationWidth > 0:
                # Erode and dilate ``highPixels`` to exclude noisy pixels.
                highPixels = ndimage.binary_opening(highPixels, structure=filterStructure)
            image[highPixels] = highThreshold[highPixels]
        if lowThreshold is not None:
            lowPixels = image < lowThreshold
            if regularizationWidth > 0:
                # Erode and dilate ``lowPixels`` to exclude noisy pixels.
                lowPixels = ndimage.binary_opening(lowPixels, structure=filterStructure)
            image[lowPixels] = lowThreshold[lowPixels]


def applyDcr(image, dcr, useInverse=False, splitSubfilters=False, splitThreshold=0.,
             doPrefilter=True, order=3):
    """Shift an image along the X and Y directions.

    Parameters
    ----------
    image : `numpy.ndarray`
        The input image to shift.
    dcr : `tuple`
        Shift calculated with ``calculateDcr``.
        Uses numpy axes ordering (Y, X).
        If ``splitSubfilters`` is set, each element is itself a `tuple`
        of two `float`, corresponding to the DCR shift at the two wavelengths.
        Otherwise, each element is a `float` corresponding to the DCR shift at
        the effective wavelength of the subfilter.
    useInverse : `bool`, optional
        Apply the shift in the opposite direction. Default: False
    splitSubfilters : `bool`, optional
        Calculate DCR for two evenly-spaced wavelengths in each subfilter,
        instead of at the midpoint. Default: False
    splitThreshold : `float`, optional
        Minimum DCR difference within a subfilter required to use
        ``splitSubfilters``
    doPrefilter : `bool`, optional
        Spline filter the image before shifting, if set. Filtering is required,
        so only set to False if the image is already filtered.
        Filtering takes ~20% of the time of shifting, so if `applyDcr` will be
        called repeatedly on the same image it is more efficient to
        precalculate the filter.
    order : `int`, optional
        The order of the spline interpolation, default is 3.

    Returns
    -------
    shiftedImage : `numpy.ndarray`
        A copy of the input image with the specified shift applied.
    """
    if doPrefilter and order > 1:
        prefilteredImage = ndimage.spline_filter(image, order=order)
    else:
        prefilteredImage = image
    if splitSubfilters:
        shiftAmp = np.max(np.abs([_dcr0 - _dcr1 for _dcr0, _dcr1 in zip(dcr[0], dcr[1])]))
        if shiftAmp >= splitThreshold:
            if useInverse:
                shift = [-1.*s for s in dcr[0]]
                shift1 = [-1.*s for s in dcr[1]]
            else:
                shift = dcr[0]
                shift1 = dcr[1]
            shiftedImage = ndimage.shift(prefilteredImage, shift, prefilter=False, order=order)
            shiftedImage += ndimage.shift(prefilteredImage, shift1, prefilter=False, order=order)
            shiftedImage /= 2.
            return shiftedImage
        else:
            # If the difference in the DCR shifts is less than the threshold,
            # then just use the average shift for efficiency.
            dcr = (np.mean(dcr[0]), np.mean(dcr[1]))
    if useInverse:
        shift = [-1.*s for s in dcr]
    else:
        shift = dcr
    shiftedImage = ndimage.shift(prefilteredImage, shift, prefilter=False, order=order)
    return shiftedImage


def calculateDcr(visitInfo, wcs, effectiveWavelength, bandwidth, dcrNumSubfilters, splitSubfilters=False,
                 bbox=None):
    """Calculate the shift in pixels of an exposure due to DCR.

    Parameters
    ----------
    visitInfo : `lsst.afw.image.VisitInfo`
        Metadata for the exposure.
    wcs : `lsst.afw.geom.SkyWcs`
        Coordinate system definition (wcs) for the exposure.
    effectiveWavelength : `float`
        The effective wavelengths of the current filter, in nanometers.
    bandwidth : `float`
        The bandwidth of the current filter, in nanometers.
    dcrNumSubfilters : `int`
        Number of sub-filters used to model chromatic effects within a band.
    splitSubfilters : `bool`, optional
        Calculate DCR for two evenly-spaced wavelengths in each subfilter,
        instead of at the midpoint. Default: False
    bbox : `lsst.afw.geom.Box2I`, optional
        Bounding box for the region of interest for evaluating the local
        pixelScale (defaults to the Sky Origin of the ``wcs`` provided if
        ``bbox`` is None).

    Returns
    -------
    dcrShift : `tuple` of two `float`
        The 2D shift due to DCR, in pixels.
        Uses numpy axes ordering (Y, X).
    """
    rotation = calculateImageParallacticAngle(visitInfo, wcs)
    dcrShift = []
    weight = [0.75, 0.25]
    for wl0, wl1 in wavelengthGenerator(effectiveWavelength, bandwidth, dcrNumSubfilters):
        # Note that diffRefractAmp can be negative, since it's relative to the
        # midpoint of the full band
        diffRefractAmp0 = differentialRefraction(wavelength=wl0, wavelengthRef=effectiveWavelength,
                                                 elevation=visitInfo.getBoresightAzAlt().getLatitude(),
                                                 observatory=visitInfo.getObservatory(),
                                                 weather=visitInfo.getWeather())
        diffRefractAmp1 = differentialRefraction(wavelength=wl1, wavelengthRef=effectiveWavelength,
                                                 elevation=visitInfo.getBoresightAzAlt().getLatitude(),
                                                 observatory=visitInfo.getObservatory(),
                                                 weather=visitInfo.getWeather())
        if bbox is not None:
            pixelScale = wcs.getPixelScale(bbox.getCenter()).asArcseconds()
        else:
            pixelScale = wcs.getPixelScale().asArcseconds()

        if splitSubfilters:
            diffRefractPix0 = diffRefractAmp0.asArcseconds()/pixelScale
            diffRefractPix1 = diffRefractAmp1.asArcseconds()/pixelScale
            diffRefractArr = [diffRefractPix0*weight[0] + diffRefractPix1*weight[1],
                              diffRefractPix0*weight[1] + diffRefractPix1*weight[0]]
            shiftX = [diffRefractPix*np.sin(rotation.asRadians()) for diffRefractPix in diffRefractArr]
            shiftY = [diffRefractPix*np.cos(rotation.asRadians()) for diffRefractPix in diffRefractArr]
            dcrShift.append(((shiftY[0], shiftX[0]), (shiftY[1], shiftX[1])))
        else:
            diffRefractAmp = (diffRefractAmp0 + diffRefractAmp1)/2.
            diffRefractPix = diffRefractAmp.asArcseconds()/pixelScale
            shiftX = diffRefractPix*np.sin(rotation.asRadians())
            shiftY = diffRefractPix*np.cos(rotation.asRadians())
            dcrShift.append((shiftY, shiftX))
    return dcrShift


def calculateImageParallacticAngle(visitInfo, wcs):
    """Calculate the total sky rotation angle of an exposure.

    Parameters
    ----------
    visitInfo : `lsst.afw.image.VisitInfo`
        Metadata for the exposure.
    wcs : `lsst.afw.geom.SkyWcs`
        Coordinate system definition (wcs) for the exposure.

    Returns
    -------
    `lsst.geom.Angle`
        The rotation of the image axis, East from North.
        Equal to the parallactic angle plus any additional rotation of the
        coordinate system.
        A rotation angle of 0 degrees is defined with
        North along the +y axis and East along the +x axis.
        A rotation angle of 90 degrees is defined with
        North along the +x axis and East along the -y axis.
    """
    parAngle = visitInfo.getBoresightParAngle().asRadians()
    cd = wcs.getCdMatrix()
    if wcs.isFlipped:
        cdAngle = (np.arctan2(-cd[0, 1], cd[0, 0]) + np.arctan2(cd[1, 0], cd[1, 1]))/2.
        rotAngle = (cdAngle + parAngle)*geom.radians
    else:
        cdAngle = (np.arctan2(cd[0, 1], -cd[0, 0]) + np.arctan2(cd[1, 0], cd[1, 1]))/2.
        rotAngle = (cdAngle - parAngle)*geom.radians
    return rotAngle


def wavelengthGenerator(effectiveWavelength, bandwidth, dcrNumSubfilters):
    """Iterate over the wavelength endpoints of subfilters.

    Parameters
    ----------
    effectiveWavelength : `float`
        The effective wavelength of the current filter, in nanometers.
    bandwidth : `float`
        The bandwidth of the current filter, in nanometers.
    dcrNumSubfilters : `int`
        Number of sub-filters used to model chromatic effects within a band.

    Yields
    ------
    `tuple` of two `float`
        The next set of wavelength endpoints for a subfilter, in nanometers.
    """
    lambdaMin = effectiveWavelength - bandwidth/2
    lambdaMax = effectiveWavelength + bandwidth/2
    wlStep = bandwidth/dcrNumSubfilters
    for wl in np.linspace(lambdaMin, lambdaMax, dcrNumSubfilters, endpoint=False):
        yield (wl, wl + wlStep)
