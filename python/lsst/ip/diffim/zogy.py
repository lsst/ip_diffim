#
# LSST Data Management System
# Copyright 2016 AURA/LSST.
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

import numpy as np

from lsst.geom import Box2I, Point2I, Extent2I
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.afw.image import ImageOrigin
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

from .imagePsfMatch import (ImagePsfMatchTask, ImagePsfMatchConfig,
                            subtractAlgorithmRegistry)

__all__ = ["ZogyTask", "ZogyConfig",
           "ZogyImagePsfMatchConfig", "ZogyImagePsfMatchTask"]


"""Tasks for performing the "Proper image subtraction" algorithm of
Zackay, et al. (2016), hereafter simply referred to as 'ZOGY (2016)'.

`ZogyTask` contains methods to perform the basic estimation of the
ZOGY diffim ``D``, its updated PSF, and the variance-normalized
likelihood image ``S_corr``. We have implemented ZOGY using the
proscribed methodology, computing all convolutions in Fourier space,
and also variants in which the convolutions are performed in real
(image) space. The former is faster and results in fewer artifacts
when the PSFs are noisy (i.e., measured, for example, via
`PsfEx`). The latter is presumed to be preferred as it can account for
masks correctly with fewer "ringing" artifacts from edge effects or
saturated stars, but noisy PSFs result in their own smaller
artifacts. Removal of these artifacts is a subject of continuing
research. Currently, we "pad" the PSFs when performing the
subtractions in real space, which reduces, but does not entirely
eliminate these artifacts.

All methods in `ZogyTask` assume template and science images are
already accurately photometrically and astrometrically registered.

`ZogyMapper` is a wrapper which runs `ZogyTask` in the
`ImageMapReduce` framework, computing of ZOGY diffim's on small,
overlapping sub-images, thereby enabling complete ZOGY diffim's which
account for spatially-varying noise and PSFs across the two input
exposures. An example of the use of this task is in the `testZogy.py`
unit test.
"""


class ZogyConfig(pexConfig.Config):
    """Configuration parameters for the ZogyTask
    """

    templateFluxScaling = pexConfig.Field(
        dtype=float,
        default=1.,
        doc="Template flux scaling factor (Fr in ZOGY paper)"
    )

    scienceFluxScaling = pexConfig.Field(
        dtype=float,
        default=1.,
        doc="Science flux scaling factor (Fn in ZOGY paper)"
    )

    scaleByCalibration = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Compute the flux normalization scaling based on the image calibration."
        "This overrides 'templateFluxScaling' and 'scienceFluxScaling'."
    )

    correctBackground = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Subtract exposure background mean to have zero expectation value."
    )

    ignoreMaskPlanes = pexConfig.ListField(
        dtype=str,
        default=("INTRP", "EDGE", "DETECTED", "SAT", "CR", "BAD", "NO_DATA", "DETECTED_NEGATIVE"),
        doc="Mask planes to ignore for statistics"
    )
    maxPsfCentroidDist = pexConfig.Field(
        dtype=float,
        default=0.2,
        doc="Maximum centroid difference allowed between the two exposure PSFs (pixels)."
    )
    doSpatialGrid = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Split the exposure and perform matching with the spatially varying PSF."
    )
    gridInnerSize = pexConfig.Field(
        dtype=float,
        default=8,
        doc="Approximate useful inner size of the grid cells in units of the "
        "estimated matching kernel size (doSpatialGrid=True only)."
    )


class ZogyTask(pipeBase.Task):
    """Task to perform ZOGY proper image subtraction. See module-level documentation for
    additional details.

    """
    ConfigClass = ZogyConfig
    _DefaultName = "imageDifferenceZogy"

    def _computeVarianceMean(self, exposure):
        """Compute the sigma-clipped mean of the variance image of ``exposure``.
        """
        statObj = afwMath.makeStatistics(exposure.variance,
                                         exposure.mask,
                                         afwMath.MEANCLIP, self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    @staticmethod
    def padCenterOriginArray(A, newShape, useInverse=False, dtype=None):
        """Zero pad an image where the origin is at the center and replace the
        origin to the corner as required by the periodic input of FFT.

        Implement also the inverse operation, crop the padding and re-center data.

        Parameters
        ----------
        A : `numpy.ndarray`
            An array to copy from.
        newShape : `tuple` of `int`
            The dimensions of the resulting array. For padding, the resulting array
            must be larger than A in each dimension. For the inverse operation this
            must be the original, before padding dimensions of the array.
        useInverse : bool, optional
            Selector of forward, add padding, operation (False)
            or its inverse, crop padding, operation (True).
        dtype: `numpy.dtype`, optional
            Dtype of output array. Values must be implicitly castable to this type.
            Use to get expected result type, e.g. single float (nympy.float32).
            If not specified, dtype is inherited from ``A``.

        Returns
        -------
        R : `numpy.ndarray`
            The padded or unpadded array with shape of `newShape` and dtype of ``dtype``.

        Notes
        -----
        For odd dimensions, the splitting is rounded to
        put the center pixel into the new corner origin (0,0). This is to be consistent
        e.g. for a dirac delta kernel that is originally located at the center pixel.


        Raises
        ------
        ValueError : ``newShape`` dimensions must be greater than or equal to the
            dimensions of ``A`` for the forward operation and less than or equal to
            for the inverse operation.
        """

        # The forward and inverse operations should round odd dimension halves at the opposite
        # sides to get the pixels back to their original positions.
        if not useInverse:
            # Forward operation: First and second halves with respect to the axes of A.
            firstHalves = [x//2 for x in A.shape]
            secondHalves = [x-y for x, y in zip(A.shape, firstHalves)]
            for d1, d2 in zip(newShape, A.shape):
                if d1 < d2:
                    raise ValueError("Newshape dimensions must be greater or equal")
        else:
            # Inverse operation: Opposite rounding
            secondHalves = [x//2 for x in newShape]
            firstHalves = [x-y for x, y in zip(newShape, secondHalves)]
            for d1, d2 in zip(newShape, A.shape):
                if d1 > d2:
                    raise ValueError("Newshape dimensions must be smaller or equal")

        if dtype is None:
            dtype = A.dtype

        R = np.zeros(newShape, dtype=dtype)
        R[-firstHalves[0]:, -firstHalves[1]:] = A[:firstHalves[0], :firstHalves[1]]
        R[:secondHalves[0], -firstHalves[1]:] = A[-secondHalves[0]:, :firstHalves[1]]
        R[:secondHalves[0], :secondHalves[1]] = A[-secondHalves[0]:, -secondHalves[1]:]
        R[-firstHalves[0]:, :secondHalves[1]] = A[:firstHalves[0], -secondHalves[1]:]
        return R

    def initializeSubImage(self, fullExp, innerBox, outerBox, noiseMeanVar, useNoise=True):
        """Initializes a sub image.

        Parameters
        ----------
        fullExp : `lsst.afw.image.Exposure`
            The full exposure to cut sub image from.
        innerBox : `lsst.geom.Box2I`
            The useful area of the calculation up to the whole bounding box of
            ``fullExp``. ``fullExp`` must contain this box.
        outerBox : `lsst.geom.Box2I`
            The overall cutting area. ``outerBox`` must be at least 1 pixel larger
            than ``inneBox`` in all directions and may not be fully contained by
            ``fullExp``.
        noiseMeanVar : `float` > 0.
            The noise variance level to initialize variance plane and to generate
            white noise for the non-overlapping region.
        useNoise : `bool`, optional
            If True, generate white noise for non-overlapping region. Otherwise,
            zero padding will be used in the non-overlapping region.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
          - ``subImg``, ``subVarImg`` : `lsst.afw.image.ImageD`
            The new sub image and its sub variance plane.

        Notes
        -----
        ``innerBox``, ``outerBox`` must be in the PARENT system of ``fullExp``.

        Supports the non-grid option when ``innerBox`` equals to the
        bounding box of ``fullExp``.
        """
        fullBox = fullExp.getBBox()
        subImg = afwImage.ImageD(outerBox, 0)
        subVarImg = afwImage.ImageD(outerBox, noiseMeanVar)
        borderBoxes = self.splitBorder(innerBox, outerBox)
        # Initialize the border region that are not fully within the exposure
        if useNoise:
            rng = np.random.Generator(
                np.random.PCG64(seed=np.array([noiseMeanVar]).view(int)))
            noiseSig = np.sqrt(noiseMeanVar)
            for box in borderBoxes:
                if not fullBox.contains(box):
                    R = subImg[box].array
                    R[...] = rng.normal(scale=noiseSig, size=R.shape)
        # Copy data to the fully contained inner region, allowing type conversion
        subImg[innerBox].array[...] = fullExp.image[innerBox].array
        subVarImg[innerBox].array[...] = fullExp.variance[innerBox].array
        # Copy data to border regions that have at least a partial overlap
        for box in borderBoxes:
            overlapBox = box.clippedTo(fullBox)
            if not overlapBox.isEmpty():
                subImg[overlapBox].array[...] = fullExp.image[overlapBox].array
                subVarImg[overlapBox].array[...] = fullExp.variance[overlapBox].array
        return pipeBase.Struct(image=subImg, variance=subVarImg)

    @staticmethod
    def estimateMatchingKernelSize(psf1, psf2):
        """Estimate the image space size of the matching kernels.

        Return ten times the larger Gaussian sigma estimate but at least
        the largest of the original psf dimensions.

        Parameters
        ----------
        psf1, psf2 : `lsst.afw.detection.Psf`
            The PSFs of the two input exposures.

        Returns
        -------
        size : `int`
            Conservative estimate for matching kernel size in pixels.
            This is the minimum padding around the inner region at each side.

        Notes
        -----
        """
        sig1 = psf1.computeShape(psf1.getAveragePosition()).getDeterminantRadius()
        sig2 = psf2.computeShape(psf2.getAveragePosition()).getDeterminantRadius()
        sig = max(sig1, sig2)
        psfBBox1 = psf1.computeBBox(psf1.getAveragePosition())
        psfBBox2 = psf2.computeBBox(psf2.getAveragePosition())
        return max(10 * sig, psfBBox1.getWidth(), psfBBox1.getHeight(),
                   psfBBox2.getWidth(), psfBBox2.getHeight())

    @staticmethod
    def splitBorder(innerBox, outerBox):
        """Split the border area around the inner box into 8 disjunct boxes.

        Parameters
        ----------
        innerBox : `lsst.geom.Box2I`
            The inner box.
        outerBox : `lsst.geom.Box2I`
            The outer box. It must be at least 1 pixel larger in each direction than the inner box.

        Returns
        -------
        resultBoxes : `list` of 8 boxes covering the edge around innerBox

        Notes
        -----
        The border boxes do not overlap. The border is covered counter clockwise
        starting from lower left corner.

        Raises
        ------
        ValueError : If ``outerBox`` is not larger than ``innerBox``.
        """
        innerBox = innerBox.dilatedBy(1)
        if not outerBox.contains(innerBox):
            raise ValueError("OuterBox must be larger by at least 1 pixel in all directions")

        # ccw sequence of corners
        o1, o2, o3, o4 = outerBox.getCorners()
        i1, i2, i3, i4 = innerBox.getCorners()
        p1 = Point2I(outerBox.minX, innerBox.minY)
        p2 = Point2I(innerBox.maxX, outerBox.minY)
        p3 = Point2I(outerBox.maxX, innerBox.maxY)
        p4 = Point2I(innerBox.minX, outerBox.maxY)

        # The 8 border boxes ccw starting from lower left
        pointPairs = ((o1, i1), (i1 + Extent2I(1, 0), p2 + Extent2I(-1, 0)), (o2, i2),
                      (i2 + Extent2I(0, 1), p3 + Extent2I(0, -1)), (o3, i3),
                      (i3 + Extent2I(-1, 0), p4 + Extent2I(1, 0)), (o4, i4),
                      (i4 + Extent2I(0, -1), p1 + Extent2I(0, 1)))
        return [Box2I(x, y, invert=True) for (x, y) in pointPairs]

    @staticmethod
    def generateGrid(imageBox, minEdgeDims, innerBoxDims, minTotalDims=None, powerOfTwo=False):
        """Generate a splitting grid for an image.

        The inner boxes cover the input image without overlap, the edges around the inner boxes do overlap
        and go beyond the image at the image edges.

        Parameters
        ----------
        imageBox : `lsst.geom.Box2I`
            Bounding box of the exposure to split.
        minEdgeDims : `lsst.geom.Extent2I`
            Minimum edge width in (x,y) directions each side.
        innerBoxDims : `lsst.geom.Extent2I`
            Minimum requested inner box dimensions (x,y).
            The actual dimensions can be larger due to rounding.
        minTotalDims: `lsst.geom.Extent2I`, optional
            If provided, minimum total outer dimensions (x,y). The edge will be increased until satisfied.
        powerOfTwo : `bool`, optional
            If True, the outer box dimensions should be rounded up to a power of 2
            by increasing the border size. This is up to 8192, above this size,
            rounding up is disabled.

        Notes
        -----
        Inner box dimensions are chosen to be as uniform as they can, remainder pixels at the edge of the
        input will be appended to the last column/row boxes.

        See diffimTests/tickets/DM-28928_spatial_grid notebooks for demonstration of this code.

        This method can be used for both PARENT and LOCAL bounding boxes.

        The outerBox dimensions are always even.

        Returns
        -------
        boxList : `list` of `lsst.pipe.base.Struct`
          ``innerBox``, ``outerBox`` : `lsst.geom.Box2I`, inner boxes and overlapping border around them.

        """
        powersOf2 = np.array([16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192])
        doubleEdgeDims = minEdgeDims * 2
        width, height = imageBox.getDimensions()
        nX = width // innerBoxDims.x  # Round down
        if nX > 0:
            innerWidth = width // nX  # Round down
        else:
            innerWidth = width
            nX = 1
        xCorners = np.zeros(nX + 1)
        xCorners[:-1] = np.arange(nX)*innerWidth + imageBox.minX
        xCorners[-1] = imageBox.endX

        nY = height // innerBoxDims.y  # Round down
        if nY > 0:
            innerHeight = height // nY  # Round down
        else:
            innerHeight = height
            nY = 1
        yCorners = np.zeros(nY + 1)
        yCorners[:-1] = np.arange(nY)*innerHeight + imageBox.minY
        yCorners[-1] = imageBox.endY

        boxes = []

        for i_y in range(nY):
            for i_x in range(nX):
                innerBox = Box2I(Point2I(xCorners[i_x], yCorners[i_y]),
                                 Point2I(xCorners[i_x + 1] - 1, yCorners[i_y + 1] - 1))

                paddedWidth = innerBox.width + doubleEdgeDims.x
                if minTotalDims is not None and paddedWidth < minTotalDims.width:
                    paddedWidth = minTotalDims.width
                if powerOfTwo:
                    i2x = np.searchsorted(powersOf2, paddedWidth, side='left')
                    if i2x < len(powersOf2):
                        paddedWidth = powersOf2[i2x]
                if paddedWidth % 2 == 1:
                    paddedWidth += 1  # Ensure total width is even

                totalXedge = paddedWidth - innerBox.width

                paddedHeight = innerBox.height + doubleEdgeDims.y
                if minTotalDims is not None and paddedHeight < minTotalDims.height:
                    paddedHeight = minTotalDims.height
                if powerOfTwo:
                    i2y = np.searchsorted(powersOf2, paddedHeight, side='left')
                    if i2y < len(powersOf2):
                        paddedHeight = powersOf2[i2y]
                if paddedHeight % 2 == 1:
                    paddedHeight += 1  # Ensure total height is even
                totalYedge = paddedHeight - innerBox.height
                outerBox = Box2I(Point2I(innerBox.minX - totalXedge//2, innerBox.minY - totalYedge//2),
                                 Extent2I(paddedWidth, paddedHeight))
                boxes.append(pipeBase.Struct(innerBox=innerBox, outerBox=outerBox))
        return boxes

    def makeSpatialPsf(self, gridPsfs):
        """Construct a CoaddPsf based on PSFs from individual sub image solutions.

        Parameters
        ----------
        gridPsfs : iterable of `lsst.pipe.base.Struct`
            Iterable of bounding boxes (``bbox``) and Psf solutions (``psf``).

        Returns
        -------
        psf : `lsst.meas.algorithms.CoaddPsf`
            A psf constructed from the PSFs of the individual subExposures.
        """
        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("weight", type="D", doc="Coadd weight")
        mycatalog = afwTable.ExposureCatalog(schema)

        # We're just using the exposure's WCS (assuming that the subExposures'
        # WCSs are the same, which they better be!).
        wcsref = self.fullExp1.getWcs()
        for i, res in enumerate(gridPsfs):
            record = mycatalog.getTable().makeRecord()
            record.setPsf(res.psf)
            record.setWcs(wcsref)
            record.setBBox(res.bbox)
            record['weight'] = 1.0
            record['id'] = i
            mycatalog.append(record)

        # create the CoaddPsf
        psf = measAlg.CoaddPsf(mycatalog, wcsref, 'weight')
        return psf

    def padAndFftImage(self, imgArr):
        """Prepare and forward FFT an image array.

        Parameters
        ----------
        imgArr : `numpy.ndarray` of `float`
            Original array. In-place modified as `numpy.nan` and `numpy.inf` are replaced by
            array mean.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            - ``imFft`` : `numpy.ndarray` of `numpy.complex`.
                FFT of image.
            - ``filtInf``, ``filtNaN`` : `numpy.ndarray` of `bool`

        Notes
        -----
        Save location of non-finite values for restoration, and replace them
        with image mean values. Re-center and zero pad array by `padCenterOriginArray`.
        """
        filtInf = np.isinf(imgArr)
        filtNaN = np.isnan(imgArr)
        imgArr[filtInf] = np.nan
        imgArr[filtInf | filtNaN] = np.nanmean(imgArr)
        self.log.debug("Replacing %s Inf and %s NaN values.",
                       np.sum(filtInf), np.sum(filtNaN))
        imgArr = self.padCenterOriginArray(imgArr, self.freqSpaceShape)
        imgArr = np.fft.fft2(imgArr)
        return pipeBase.Struct(imFft=imgArr, filtInf=filtInf, filtNaN=filtNaN)

    def removeNonFinitePixels(self, imgArr):
        """Replace non-finite pixel values in-place.

        Save the locations of non-finite values for restoration, and replace them
        with image mean values.

        Parameters
        ----------
        imgArr : `numpy.ndarray` of `float`
            The image array. Non-finite values are replaced in-place in this array.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            - ``filtInf``, ``filtNaN`` : `numpy.ndarray` of `bool`
                The filter of the pixel values that were inf or nan.
        """
        filtInf = np.isinf(imgArr)
        filtNaN = np.isnan(imgArr)
        # Masked edge and bad pixels could also be removed here in the same way
        # in the future
        imgArr[filtInf] = np.nan
        imgArr[filtInf | filtNaN] = np.nanmean(imgArr)
        self.log.debug("Replacing %s Inf and %s NaN values.",
                       np.sum(filtInf), np.sum(filtNaN))
        return pipeBase.Struct(filtInf=filtInf, filtNaN=filtNaN)

    def inverseFftAndCropImage(self, imgArr, origSize, filtInf=None, filtNaN=None, dtype=None):
        """Inverse FFT and crop padding from image array.

        Parameters
        ----------
        imgArr : `numpy.ndarray` of `numpy.complex`
            Fourier space array representing a real image.

        origSize : `tuple` of `int`
            Original unpadded shape tuple of the image to be cropped to.

        filtInf, filtNan : `numpy.ndarray` of bool or int, optional
            If specified, they are used as index arrays for ``result`` to set values to
            `numpy.inf` and `numpy.nan` respectively at these positions.

        dtype : `numpy.dtype`, optional
            Dtype of result array to cast return values to implicitly. This is to
            spare one array copy operation at reducing double precision to single.
            If `None` result inherits dtype of `imgArr`.

        Returns
        -------
        result : `numpy.ndarray` of `dtype`
        """
        imgNew = np.fft.ifft2(imgArr)
        imgNew = imgNew.real
        imgNew = self.padCenterOriginArray(imgNew, origSize, useInverse=True, dtype=dtype)
        if filtInf is not None:
            imgNew[filtInf] = np.inf
        if filtNaN is not None:
            imgNew[filtNaN] = np.nan
        return imgNew

    @staticmethod
    def computePsfAtCenter(exposure):
        """Computes the PSF image at the bbox center point.

        This may be at a fractional pixel position.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure with psf.

        Returns
        -------
        psfImg : `lsst.afw.image.Image`
            Calculated psf image.
        """
        pBox = exposure.getBBox()
        cen = pBox.getCenter()
        psf = exposure.getPsf()
        psfImg = psf.computeKernelImage(cen)  # Centered and normed
        return psfImg

    @staticmethod
    def subtractImageMean(image, mask, statsControl):
        """In-place subtraction of sigma-clipped mean of the image.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to manipulate. Its sigma clipped mean is in-place subtracted.

        mask : `lsst.afw.image.Mask`
            Mask to use for ignoring pixels.

        statsControl : `lsst.afw.math.StatisticsControl`
            Config of sigma clipped mean statistics calculation.

        Returns
        -------
        None

        Raises
        ------
        ValueError : If image mean is nan.
        """
        statObj = afwMath.makeStatistics(image, mask,
                                         afwMath.MEANCLIP, statsControl)
        mean = statObj.getValue(afwMath.MEANCLIP)
        if not np.isnan(mean):
            image -= mean
        else:
            raise ValueError("Image mean is NaN.")

    def prepareFullExposure(self, exposure1, exposure2, correctBackground=False):
        """Performs calculations that apply to the full exposures once only.

        Parameters
        ----------

        exposure1, exposure2 : `lsst.afw.image.Exposure`
            The input exposures. Copies are made for internal calculations.

        correctBackground : `bool`, optional
            If True, subtracts sigma-clipped mean of exposures. The algorithm
            assumes zero expectation value at background pixels.

        Returns
        -------
        None

        Notes
        -----
        Set a number of instance fields with pre-calculated values.

        Raises
        ------
        ValueError : If photometric calibrations are not available while
            ``config.scaleByCalibration`` equals True.
        """
        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.Mask.getPlaneBitMask(
            self.config.ignoreMaskPlanes))

        exposure1 = exposure1.clone()
        exposure2 = exposure2.clone()
        # Fallback values if sub exposure variance calculation is problematic
        sig1 = np.sqrt(self._computeVarianceMean(exposure1))
        self.fullExpVar1 = sig1*sig1
        sig2 = np.sqrt(self._computeVarianceMean(exposure2))
        self.fullExpVar2 = sig2*sig2

        # If 'scaleByCalibration' is True then these norms are overwritten
        if self.config.scaleByCalibration:
            calibObj1 = exposure1.getPhotoCalib()
            calibObj2 = exposure2.getPhotoCalib()
            if calibObj1 is None or calibObj2 is None:
                raise ValueError("Photometric calibrations are not available for both exposures.")
            mImg1 = calibObj1.calibrateImage(exposure1.maskedImage)
            mImg2 = calibObj2.calibrateImage(exposure2.maskedImage)
            self.F1 = 1.
            self.F2 = 1.
        else:
            self.F1 = self.config.templateFluxScaling  # default is 1
            self.F2 = self.config.scienceFluxScaling  # default is 1
            mImg1 = exposure1.maskedImage
            mImg2 = exposure2.maskedImage

        # mImgs can be in-place modified
        if correctBackground:
            self.subtractImageMean(mImg1.image, mImg1.mask, self.statsControl)
            self.subtractImageMean(mImg2.image, mImg2.mask, self.statsControl)

        # Determine border size
        self.borderSize = self.estimateMatchingKernelSize(exposure1.getPsf(), exposure2.getPsf())
        self.log.debug("Minimum padding border size: %s pixels", self.borderSize)
        # Remove non-finite values from the images in-place
        self.filtsImg1 = self.removeNonFinitePixels(mImg1.image.array)
        self.filtsImg2 = self.removeNonFinitePixels(mImg2.image.array)
        self.filtsVar1 = self.removeNonFinitePixels(mImg1.variance.array)
        self.filtsVar2 = self.removeNonFinitePixels(mImg2.variance.array)

        exposure1.maskedImage = mImg1
        exposure2.maskedImage = mImg2

        self.fullExp1 = exposure1
        self.fullExp2 = exposure2

    def prepareSubExposure(self, localCutout, psf1=None, psf2=None, sig1=None, sig2=None):
        """Perform per-sub exposure preparations.

        Parameters
        ----------
        sig1, sig2 : `float`, optional
            For debug purposes only, copnsider that the image
            may already be rescaled by the photometric calibration.
        localCutout : `lsst.pipe.base.Struct`
            - innerBox, outerBox: `lsst.geom.Box2I` LOCAL inner and outer boxes
        psf1, psf2 : `lsst.afw.detection.Psf`, optional
            If specified, use given psf as the sub exposure psf. For debug purposes.
        sig1, sig2 : `float`, optional
            If specified, use value as the sub-exposures' background noise sigma value.

        Returns
        -------
        None

        """
        self.log.debug("Processing LOCAL cell w/ inner box:%s, outer box:%s",
                       localCutout.innerBox, localCutout.outerBox)
        # The PARENT origin cutout boxes for the two exposures
        self.cutBoxes1 = pipeBase.Struct(
            innerBox=localCutout.innerBox.shiftedBy(Extent2I(self.fullExp1.getXY0())),
            outerBox=localCutout.outerBox.shiftedBy(Extent2I(self.fullExp1.getXY0())))
        self.cutBoxes2 = pipeBase.Struct(
            innerBox=localCutout.innerBox.shiftedBy(Extent2I(self.fullExp2.getXY0())),
            outerBox=localCutout.outerBox.shiftedBy(Extent2I(self.fullExp2.getXY0())))
        # The sub-exposure views of the useful inner area of this grid cell
        innerSubExp1 = self.fullExp1[self.cutBoxes1.innerBox]
        innerSubExp2 = self.fullExp2[self.cutBoxes2.innerBox]
        if psf1 is None:
            self.subExpPsf1 = self.computePsfAtCenter(innerSubExp1)
        else:
            self.subExpPsf1 = psf1
        if psf2 is None:
            self.subExpPsf2 = self.computePsfAtCenter(innerSubExp2)
        else:
            self.subExpPsf2 = psf2
        self.checkCentroids(self.subExpPsf1.array, self.subExpPsf2.array)
        psfBBox1 = self.subExpPsf1.getBBox()
        psfBBox2 = self.subExpPsf2.getBBox()
        self.psfShape1 = (psfBBox1.getHeight(), psfBBox1.getWidth())
        self.psfShape2 = (psfBBox2.getHeight(), psfBBox2.getWidth())
        # sig1 and sig2  should not be set externally, just for debug purpose
        if sig1 is None:
            sig1 = np.sqrt(self._computeVarianceMean(innerSubExp1))
        if sig1 > 0.:  # Not zero and not nan
            self.subExpVar1 = sig1*sig1
        else:
            self.subExpVar1 = self.fullExpVar1
        if sig2 is None:
            sig2 = np.sqrt(self._computeVarianceMean(innerSubExp2))
        if sig2 > 0.:  # Not zero and not nan
            self.subExpVar2 = sig2*sig2
        else:
            self.subExpVar2 = self.fullExpVar2
        self.freqSpaceShape = (localCutout.outerBox.getHeight(), localCutout.outerBox.getWidth())

        self.subImg1 = self.initializeSubImage(
            self.fullExp1, self.cutBoxes1.innerBox, self.cutBoxes1.outerBox,
            self.subExpVar1, useNoise=True)
        self.subImg2 = self.initializeSubImage(
            self.fullExp2, self.cutBoxes2.innerBox, self.cutBoxes2.outerBox,
            self.subExpVar2, useNoise=True)

        D = self.padCenterOriginArray(self.subImg1.image.array, self.freqSpaceShape)
        self.subImgFft1 = np.fft.fft2(D)
        D = self.padCenterOriginArray(self.subImg1.variance.array, self.freqSpaceShape)
        self.subVarImgFft1 = np.fft.fft2(D)

        D = self.padCenterOriginArray(self.subImg2.image.array, self.freqSpaceShape)
        self.subImgFft2 = np.fft.fft2(D)
        D = self.padCenterOriginArray(self.subImg2.variance.array, self.freqSpaceShape)
        self.subVarImgFft2 = np.fft.fft2(D)

        D = self.padCenterOriginArray(self.subExpPsf1.array, self.freqSpaceShape)
        self.psfFft1 = np.fft.fft2(D)
        D = self.padCenterOriginArray(self.subExpPsf2.array, self.freqSpaceShape)
        self.psfFft2 = np.fft.fft2(D)

    @staticmethod
    def pixelSpaceSquare(D):
        """Square the argument in pixel space.

        Parameters
        ----------
        D : 2D `numpy.ndarray` of `numpy.complex`
            Fourier transform of a real valued array.

        Returns
        -------
        R : `numpy.ndarray` of `numpy.complex`

        Notes
        -----
        ``D`` is to be inverse Fourier transformed, squared and then
        forward Fourier transformed again, i.e. an autoconvolution in Fourier space.
        This operation is not distributive over multiplication.
        ``pixelSpaceSquare(A*B) != pixelSpaceSquare(A)*pixelSpaceSquare(B)``
        """
        R = np.real(np.fft.ifft2(D))
        R *= R
        R = np.fft.fft2(R)
        return R

    @staticmethod
    def getCentroid(A):
        """Calculate the centroid coordinates of a 2D array.

        Parameters
        ----------
        A : 2D `numpy.ndarray` of `float`
            The input array. Must not be all exact zero.

        Notes
        -----
        Calculates the centroid as if the array represented a 2D geometrical shape with
        weights per cell, allowing for "negative" weights. If sum equals to exact (float) zero,
        calculates centroid of absolute value array.

        The geometrical center is defined as (0,0), independently of the array shape.
        For an odd dimension, this is the center of the center pixel,
        for an even dimension, this is between the two center pixels.

        Returns
        -------
        ycen, xcen : `tuple` of `float`

        """
        s = np.sum(A)
        if s == 0.:
            A = np.fabs(A)
            s = np.sum(A)
        w = np.arange(A.shape[0], dtype=float) - (A.shape[0] - 1.)/2
        ycen = np.sum(w[:, np.newaxis]*A)/s
        w = np.arange(A.shape[1], dtype=float) - (A.shape[1] - 1.)/2
        xcen = np.sum(w[np.newaxis, :]*A)/s

        return ycen, xcen

    def checkCentroids(self, psfArr1, psfArr2):
        """Check whether two PSF array centroids' distance is within tolerance.

        Parameters
        ----------
        psfArr1, psfArr2 : `numpy.ndarray` of `float`
            Input PSF arrays to check.

        Returns
        -------
        None

        Raises
        ------
        ValueError:
            Centroid distance exceeds `config.maxPsfCentroidDist` pixels.
        """
        yc1, xc1 = self.getCentroid(psfArr1)
        yc2, xc2 = self.getCentroid(psfArr2)
        dy = yc2 - yc1
        dx = xc2 - xc1
        if dy*dy + dx*dx > self.config.maxPsfCentroidDist*self.config.maxPsfCentroidDist:
            raise ValueError(
                f"PSF centroids are offset by more than {self.config.maxPsfCentroidDist:.2f} pixels.")

    def calculateFourierDiffim(self, psf1, im1, varPlane1, F1, varMean1,
                               psf2, im2, varPlane2, F2, varMean2, calculateScore=True):
        """Convolve and subtract two images in Fourier space.

        Calculate the ZOGY proper difference image, score image and their PSFs.
        All input and output arrays are in Fourier space.

        Parameters
        ----------
        psf1, psf2 : `numpy.ndarray`, (``self.freqSpaceShape``,)
            Psf arrays. Must be already in Fourier space.
        im1, im2 : `numpy.ndarray`, (``self.freqSpaceShape``,)
            Image arrays. Must be already in Fourier space.
        varPlane1, varPlane2 : `numpy.ndarray`, (``self.freqSpaceShape``,)
            Variance plane arrays respectively. Must be already in Fourier space.
        varMean1, varMean2 : `numpy.float` > 0.
            Average per-pixel noise variance in im1, im2 respectively. Used as weighing
            of input images. Must be greater than zero.
        F1, F2 : `numpy.float` > 0.
            Photometric scaling of the images. See eqs. (5)--(9)
        calculateScore : `bool`, optional
            If True (default), calculate and return the detection significance (score) image.
            Otherwise, these return fields are `None`.

        Returns
        -------
        result : `pipe.base.Struct`
            All arrays are in Fourier space and have shape ``self.freqSpaceShape``.

            ``Fd``
                Photometric level of ``D`` (`float`).
            ``D``
                The difference image (`numpy.ndarray` [`numpy.complex`]).
            ``varplaneD``
                Variance plane of ``D`` (`numpy.ndarray` [`numpy.complex`]).
            ``Pd``
                PSF of ``D`` (`numpy.ndarray` [`numpy.complex`]).
            ``S``
                Significance (score) image (`numpy.ndarray` [`numpy.complex`] or `None`).
            ``varplaneS``
                Variance plane of ``S`` ((`numpy.ndarray` [`numpy.complex`] or `None`).
            ``Ps``
                PSF of ``S`` (`numpy.ndarray` [`numpy.complex`]).

        Notes
        -----
        All array inputs and outputs are Fourier-space images with shape of
        `self.freqSpaceShape` in this method.

        ``varMean1``, ``varMean2`` quantities are part of the noise model and not to be confused
        with the variance of image frequency components or with ``varPlane1``, ``varPlane2`` that
        are the Fourier transform of the variance planes.
        """
        var1F2Sq = varMean1*F2*F2
        var2F1Sq = varMean2*F1*F1
        # We need reals for comparison, also real operations are usually faster
        psfAbsSq1 = np.real(np.conj(psf1)*psf1)
        psfAbsSq2 = np.real(np.conj(psf2)*psf2)
        FdDenom = np.sqrt(var1F2Sq + var2F1Sq)  # one number

        # Secure positive limit to avoid floating point operations resulting in exact zero
        tiny = np.finfo(psf1.dtype).tiny * 100
        sDenom = var1F2Sq*psfAbsSq2 + var2F1Sq*psfAbsSq1  # array, eq. (12)
        # Frequencies where both psfs are too close to zero.
        # We expect this only in cases when psf1, psf2 are identical,
        # and either having very well sampled Gaussian tails
        # or having "edges" such that some sinc-like zero crossings are found at symmetry points
        #
        # if sDenom < tiny then it can be == 0. -> `denom` = 0. and 0/0 occur at `c1` , `c2`
        # if we keep SDenom = tiny, denom ~ O(sqrt(tiny)), Pd ~ O(sqrt(tiny)), S ~ O(sqrt(tiny)*tiny) == 0
        # Where S = 0 then Pd = 0 and D should still yield the same variance ~ O(1)
        # For safety, we set S = 0 explicitly, too, though it should be unnecessary.
        fltZero = sDenom < tiny
        nZero = np.sum(fltZero)
        self.log.debug("There are %s frequencies where both FFTd PSFs are close to zero.", nZero)
        if nZero > 0:
            # We expect only a small fraction of such frequencies
            fltZero = np.nonzero(fltZero)  # Tuple of index arrays
            sDenom[fltZero] = tiny  # Avoid division problem but overwrite result anyway
        denom = np.sqrt(sDenom)  # array, eq. (13)

        c1 = F2*psf2/denom
        c2 = F1*psf1/denom
        if nZero > 0:
            c1[fltZero] = F2/FdDenom
            c2[fltZero] = F1/FdDenom
        D = c1*im1 - c2*im2  # Difference image eq. (13)
        varPlaneD = self.pixelSpaceSquare(c1)*varPlane1 + self.pixelSpaceSquare(c2)*varPlane2  # eq. (26)

        Pd = FdDenom*psf1*psf2/denom  # Psf of D eq. (14)
        if nZero > 0:
            Pd[fltZero] = 0

        Fd = F1*F2/FdDenom  # Flux scaling of D eq. (15)
        if calculateScore:
            c1 = F1*F2*F2*np.conj(psf1)*psfAbsSq2/sDenom
            c2 = F2*F1*F1*np.conj(psf2)*psfAbsSq1/sDenom
            if nZero > 0:
                c1[fltZero] = 0
                c2[fltZero] = 0
            S = c1*im1 - c2*im2  # eq. (12)
            varPlaneS = self.pixelSpaceSquare(c1)*varPlane1 + self.pixelSpaceSquare(c2)*varPlane2
            Ps = np.conj(Pd)*Pd  # eq. (17) Source detection expects a PSF
        else:
            S = None
            Ps = None
            varPlaneS = None
        return pipeBase.Struct(D=D, Pd=Pd, varPlaneD=varPlaneD, Fd=Fd,
                               S=S, Ps=Ps, varPlaneS=varPlaneS)

    @staticmethod
    def calculateMaskPlane(mask1, mask2, effPsf1=None, effPsf2=None):
        """Calculate the mask plane of the difference image.

        Parameters
        ----------
        mask1, maks2 : `lsst.afw.image.Mask`
            Mask planes of the two exposures.


        Returns
        -------
        diffmask : `lsst.afw.image.Mask`
            Mask plane for the subtraction result.

        Notes
        -----
        TODO DM-25174 : Specification of effPsf1, effPsf2 are not yet supported.
        """

        # mask1 x effPsf2 | mask2 x effPsf1
        if effPsf1 is not None or effPsf2 is not None:
            # TODO: DM-25174 effPsf1, effPsf2: the effective psf for cross-blurring.
            # We need a "size" approximation of the c1 and c2 coefficients to make effPsfs
            # Also convolution not yet supports mask-only operation
            raise NotImplementedError("Mask plane only 'convolution' operation is not yet supported")
        R = mask1.clone()
        R |= mask2
        return R

    @staticmethod
    def makeKernelPsfFromArray(A):
        """Create a non spatially varying PSF from a `numpy.ndarray`.

        Parameters
        ----------
        A : `numpy.ndarray`
            2D array to use as the new psf image. The pixels are copied.

        Returns
        -------
        psfNew : `lsst.meas.algorithms.KernelPsf`
            The constructed PSF.
        """
        psfImg = afwImage.ImageD(A.astype(np.float64, copy=True), deep=False)
        psfNew = measAlg.KernelPsf(afwMath.FixedKernel(psfImg))
        return psfNew

    def pasteSubDiffImg(self, ftDiff, diffExp, scoreExp=None):
        """Paste sub image results back into result Exposure objects.

        Parameters
        ----------
        ftDiff : `lsst.pipe.base.Struct`
            Result struct by `calculateFourierDiffim`.
        diffExp : `lsst.afw.image.Exposure`
            The result exposure to paste into the sub image result.
            Must be dimensions and dtype of ``self.fullExp1``.
        scoreExp : `lsst.afw.image.Exposure` or `None`
            The result score exposure to paste into the sub image result.
            Must be dimensions and dtype of ``self.fullExp1``.
            If `None`, the score image results are disregarded.

        Returns
        -------
        None

        Notes
        -----
        The PSF of the score image is just to make the score image resemble a
        regular exposure and to study the algorithm performance.

        Add an entry to the ``self.gridPsfs`` list.

        gridPsfs : `list` of `lsst.pipe.base.Struct`
            - ``bbox`` : `lsst.geom.Box2I`
                The inner region of the grid cell.
            - ``Pd`` :  `lsst.meas.algorithms.KernelPsf`
                The diffim PSF in this cell.
            - ``Ps`` :  `lsst.meas.algorithms.KernelPsf` or `None`
                The score image PSF in this cell or `None` if the score
                image was not calculated.
        """
        D = self.inverseFftAndCropImage(
            ftDiff.D, self.freqSpaceShape, dtype=self.fullExp1.image.dtype)
        varPlaneD = self.inverseFftAndCropImage(
            ftDiff.varPlaneD, self.freqSpaceShape, dtype=self.fullExp1.variance.dtype)
        Pd = self.inverseFftAndCropImage(
            ftDiff.Pd, self.psfShape1, dtype=self.subExpPsf1.dtype)
        sumPd = np.sum(Pd)
        # If this is smaller than 1. it is an indicator that it does not fit its original dimensions
        self.log.info("Pd sum before normalization: %.3f", sumPd)
        Pd /= sumPd
        # Convert Pd into a Psf object
        Pd = self.makeKernelPsfFromArray(Pd)

        xy0 = self.cutBoxes1.outerBox.getMin()
        # D is already converted back to dtype of fullExp1
        # Encapsulate D simply into an afwImage.Image for correct inner-outer box handling
        imgD = afwImage.Image(D, deep=False, xy0=xy0, dtype=self.fullExp1.image.dtype)
        diffExp.image[self.cutBoxes1.innerBox] = imgD[self.cutBoxes1.innerBox]
        imgVarPlaneD = afwImage.Image(varPlaneD, deep=False, xy0=xy0,
                                      dtype=self.fullExp1.variance.dtype)
        diffExp.variance[self.cutBoxes1.innerBox] = imgVarPlaneD[self.cutBoxes1.innerBox]
        diffExp.mask[self.cutBoxes1.innerBox] = self.calculateMaskPlane(
            self.fullExp1.mask[self.cutBoxes1.innerBox], self.fullExp2.mask[self.cutBoxes2.innerBox])

        # Calibrate the image; subimages on the grid must be on the same photometric scale
        # Now the calibration object will be 1. everywhere
        diffExp.maskedImage[self.cutBoxes1.innerBox] /= ftDiff.Fd

        if ftDiff.S is not None and scoreExp is not None:
            S = self.inverseFftAndCropImage(
                ftDiff.S, self.freqSpaceShape, dtype=self.fullExp1.image.dtype)
            varPlaneS = self.inverseFftAndCropImage(
                ftDiff.varPlaneS, self.freqSpaceShape, dtype=self.fullExp1.variance.dtype)

            imgS = afwImage.Image(S, deep=False, xy0=xy0, dtype=self.fullExp1.image.dtype)
            imgVarPlaneS = afwImage.Image(varPlaneS, deep=False, xy0=xy0,
                                          dtype=self.fullExp1.variance.dtype)
            scoreExp.image[self.cutBoxes1.innerBox] = imgS[self.cutBoxes1.innerBox]
            scoreExp.variance[self.cutBoxes1.innerBox] = imgVarPlaneS[self.cutBoxes1.innerBox]

            # PSF of S
            Ps = self.inverseFftAndCropImage(ftDiff.Ps, self.psfShape1, dtype=self.subExpPsf1.dtype)
            sumPs = np.sum(Ps)
            self.log.info("Ps sum before normalization: %.3f", sumPs)
            Ps /= sumPs

            # TODO DM-23855 : Additional score image corrections may be done here

            scoreExp.mask[self.cutBoxes1.innerBox] = diffExp.mask[self.cutBoxes1.innerBox]
            # Convert Ps into a Psf object
            Ps = self.makeKernelPsfFromArray(Ps)
        else:
            Ps = None
        self.gridPsfs.append(pipeBase.Struct(bbox=self.cutBoxes1.innerBox, Pd=Pd, Ps=Ps))

    def finishResultExposures(self, diffExp, scoreExp=None):
        """Perform final steps on the full difference exposure result.

        Set photometric calibration, psf properties of the exposures.

        Parameters
        ----------
        diffExp : `lsst.afw.image.Exposure`
            The result difference image exposure to finalize.
        scoreExp : `lsst.afw.image.Exposure` or `None`
            The result score exposure to finalize.

        Returns
        -------
        None.
        """
        # Set Calibration and PSF of the result exposures
        calibOne = afwImage.PhotoCalib(1.)
        diffExp.setPhotoCalib(calibOne)
        # Create the spatially varying PSF and set it for the diffExp
        # Set the PSF of this subExposure
        if len(self.gridPsfs) > 1:
            diffExp.setPsf(
                self.makeSpatialPsf(
                    pipeBase.Struct(bbox=x.bbox, psf=x.Pd) for x in self.gridPsfs
                ))
            if scoreExp is not None:
                scoreExp.setPsf(
                    self.makeSpatialPsf(
                        pipeBase.Struct(bbox=x.bbox, psf=x.Ps) for x in self.gridPsfs
                    ))
        else:
            # We did not have a grid, use the result psf without
            # making a CoaddPsf
            diffExp.setPsf(self.gridPsfs[0].Pd)
            if scoreExp is not None:
                scoreExp.setPsf(self.gridPsfs[0].Ps)

        # diffExp.setPsf(self.makeKernelPsfFromArray(Pd))
        if scoreExp is not None:
            scoreExp.setPhotoCalib(calibOne)
            # Zero score exposure where its variance is zero or the inputs are non-finite
            flt = (self.filtsImg1.filtInf | self.filtsImg2.filtInf
                   | self.filtsImg1.filtNaN | self.filtsImg2.filtNaN
                   | self.filtsVar1.filtInf | self.filtsVar2.filtInf
                   | self.filtsVar1.filtNaN | self.filtsVar2.filtNaN)
            # Ensure that no division by 0 occurs in S/sigma(S).
            # S is set to be always finite, 0 where pixels non-finite
            tiny = np.finfo(scoreExp.variance.dtype).tiny * 100
            flt = np.logical_or(flt, scoreExp.variance.array < tiny)
            # Don't set variance to tiny.
            # It becomes 0 in case of conversion to single precision.
            # Set variance to 1, indicating that zero is in units of "sigmas" already.
            scoreExp.variance.array[flt] = 1
            scoreExp.image.array[flt] = 0

    def run(self, exposure1, exposure2, calculateScore=True):
        """Task entry point to perform the zogy subtraction
        of ``exposure1-exposure2``.

        Parameters
        ----------
        exposure1, exposure2 : `lsst.afw.image.Exposure`
            Two exposures warped and matched into matching pixel dimensions.
        calculateScore : `bool`, optional
            If True (default), calculate the score image and return in ``scoreExp``.


        Returns
        -------
        resultName : `lsst.pipe.base.Struct`
            - ``diffExp`` : `lsst.afw.image.Exposure`
                The Zogy difference exposure (``exposure1-exposure2``).
            - ``scoreExp`` : `lsst.afw.image.Exposure` or `None`
                The Zogy significance or score (S) exposure if ``calculateScore==True``.
            - ``ftDiff`` : `lsst.pipe.base.Struct`
                Lower level return struct by `calculateFourierDiffim` with added
                fields from the task instance. For debug purposes.

        Notes
        -----

        ``diffExp`` and ``scoreExp`` always inherit their metadata from
        ``exposure1`` (e.g. dtype, bbox, wcs).

        The score image (``S``) is defined in the ZOGY paper as the detection
        statistic value at each pixel. In the ZOGY image model, the input images
        have uniform variance noises and thus ``S`` has uniform per pixel
        variance (though it is not scaled to 1). In Section 3.3 of the paper,
        there are "corrections" defined to the score image to correct the
        significance values for some deviations from the image model. The first
        of these corrections is the calculation of the *variance plane* of ``S``
        allowing for different per pixel variance values by following the
        overall convolution operation on the pixels of the input images. ``S``
        scaled (divided) by its corrected per pixel noise is referred as
        ``Scorr`` in the paper.

        In the current implementation, ``scoreExp`` contains ``S`` in its image
        plane and the calculated (non-uniform) variance plane of ``S`` in its
        variance plane. ``scoreExp`` can be used directly for source detection
        as a likelihood image by respecting its variance plane or can be divided
        by the square root of the variance plane to scale detection significance
        values into units of sigma. ``S`` should be interpreted as a detection
        likelihood directly on a per-pixel basis. The calculated PSF
        of ``S`` is merely an indication how much the input PSFs localize point
        sources.

        TODO DM-23855 : Implement further correction tags to the variance of
        ``scoreExp``. As of DM-25174 it is not determined how important these
        further correction tags are.
        """
        # We use the dimensions of the 1st image only in the code
        if exposure1.getDimensions() != exposure2.getDimensions():
            raise ValueError("Exposure dimensions do not match ({} != {} )".format(
                             exposure1.getDimensions(), exposure2.getDimensions()))

        self.prepareFullExposure(exposure1, exposure2, correctBackground=self.config.correctBackground)
        # Do not use the exposure1, exposure2 input arguments from here
        exposure1 = None
        exposure2 = None
        if self.config.doSpatialGrid:
            gridBoxes = self.generateGrid(
                self.fullExp1.getBBox(ImageOrigin.LOCAL), Extent2I(self.borderSize, self.borderSize),
                Extent2I(Extent2I(self.borderSize, self.borderSize) * self.config.gridInnerSize),
                powerOfTwo=True)
        else:
            gridBoxes = self.generateGrid(
                self.fullExp1.getBBox(ImageOrigin.LOCAL), Extent2I(self.borderSize, self.borderSize),
                self.fullExp1.getBBox().getDimensions(), powerOfTwo=True)

        diffExp = self.fullExp1.clone()
        if calculateScore:
            scoreExp = self.fullExp1.clone()
        else:
            scoreExp = None
        self.gridPsfs = []
        # Loop through grid boxes
        for boxPair in gridBoxes:
            self.prepareSubExposure(boxPair)  # Extract sub images and fft
            ftDiff = self.calculateFourierDiffim(
                self.psfFft1, self.subImgFft1, self.subVarImgFft1, self.F1, self.subExpVar1,
                self.psfFft2, self.subImgFft2, self.subVarImgFft2, self.F2, self.subExpVar2,
                calculateScore=calculateScore)
            self.pasteSubDiffImg(ftDiff, diffExp, scoreExp)  # Paste back result
        self.finishResultExposures(diffExp, scoreExp)
        # Add debug info from the task instance
        ftDiff.freqSpaceShape = self.freqSpaceShape  # The outer shape of the last grid cell
        ftDiff.psfShape1 = self.psfShape1  # The psf image shape in exposure1
        ftDiff.psfShape2 = self.psfShape2  # The psf image shape in exposure2
        ftDiff.borderSize = self.borderSize  # The requested padding around the inner region
        return pipeBase.Struct(diffExp=diffExp,
                               scoreExp=scoreExp,
                               ftDiff=ftDiff)


class ZogyImagePsfMatchConfig(ImagePsfMatchConfig):
    """Config for the ZogyImagePsfMatchTask"""

    zogyConfig = pexConfig.ConfigField(
        dtype=ZogyConfig,
        doc='ZogyTask config to use',
    )


class ZogyImagePsfMatchTask(ImagePsfMatchTask):
    """Task to perform Zogy PSF matching and image subtraction.

    This class inherits from ImagePsfMatchTask to contain the _warper
    subtask and related methods.
    """

    ConfigClass = ZogyImagePsfMatchConfig

    def __init__(self, *args, **kwargs):
        ImagePsfMatchTask.__init__(self, *args, **kwargs)

    def run(self, scienceExposure, templateExposure, doWarping=True):
        """Register, PSF-match, and subtract two Exposures, ``scienceExposure - templateExposure``
        using the ZOGY algorithm.

        Parameters
        ----------
        templateExposure : `lsst.afw.image.Exposure`
            exposure to be warped to scienceExposure.
        scienceExposure : `lsst.afw.image.Exposure`
            reference Exposure.
        doWarping : `bool`
            what to do if templateExposure's and scienceExposure's WCSs do not match:
            - if True then warp templateExposure to match scienceExposure
            - if False then raise an Exception

        Notes
        -----
        Do the following, in order:
            - Warp templateExposure to match scienceExposure, if their WCSs do not already match
            - Compute subtracted exposure ZOGY image subtraction algorithm on the two exposures

        This is the new entry point of the task as of DM-25115.

        Returns
        -------
        results : `lsst.pipe.base.Struct` containing these fields:
            - subtractedExposure: `lsst.afw.image.Exposure`
                The subtraction result.
            - warpedExposure: `lsst.afw.image.Exposure` or `None`
                templateExposure after warping to match scienceExposure
        """

        if not self._validateWcs(scienceExposure, templateExposure):
            if doWarping:
                self.log.info("Warping templateExposure to scienceExposure")
                xyTransform = afwGeom.makeWcsPairTransform(templateExposure.getWcs(),
                                                           scienceExposure.getWcs())
                psfWarped = measAlg.WarpedPsf(templateExposure.getPsf(), xyTransform)
                templateExposure = self._warper.warpExposure(
                    scienceExposure.getWcs(), templateExposure, destBBox=scienceExposure.getBBox())
                templateExposure.setPsf(psfWarped)
            else:
                raise RuntimeError("Input images are not registered. Consider setting doWarping=True.")

        config = self.config.zogyConfig
        task = ZogyTask(config=config)
        results = task.run(scienceExposure, templateExposure)
        results.warpedExposure = templateExposure
        return results

    def subtractExposures(self, templateExposure, scienceExposure, *args):
        raise NotImplementedError

    def subtractMaskedImages(self, templateExposure, scienceExposure, *args):
        raise NotImplementedError


subtractAlgorithmRegistry.register('zogy', ZogyImagePsfMatchTask)
