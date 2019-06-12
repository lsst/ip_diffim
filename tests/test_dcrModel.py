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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.

from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, Angle
from astropy.time import Time
import numpy as np
from scipy import ndimage
import unittest

from astro_metadata_translator import makeObservationInfo
from lsst.afw.coord.refraction import differentialRefraction
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.math as afwMath
from lsst.geom import arcseconds, degrees, radians, arcminutes
from lsst.ip.diffim.dcrModel import (DcrModel, calculateDcr, calculateImageParallacticAngle,
                                     applyDcr, wavelengthGenerator)
from lsst.obs.base import MakeRawVisitInfoViaObsInfo
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests


# Our calculation of hour angle and parallactic angle ignore precession
# and nutation, so calculations depending on these are not precise. DM-20133
coordinateTolerance = 1.*arcminutes


class DcrModelTestTask(lsst.utils.tests.TestCase):
    """A test case for the DCR-aware image coaddition algorithm.

    Attributes
    ----------
    bbox : `lsst.afw.geom.Box2I`
        Bounding box of the test model.
    bufferSize : `int`
        Distance from the inner edge of the bounding box
        to avoid placing test sources in the model images.
    dcrNumSubfilters : int
        Number of sub-filters used to model chromatic effects within a band.
    lambdaEff : `float`
        Effective wavelength of the full band.
    lambdaMax : `float`
        Maximum wavelength where the relative throughput
        of the band is greater than 1%.
    lambdaMin : `float`
        Minimum wavelength where the relative throughput
        of the band is greater than 1%.
    mask : `lsst.afw.image.Mask`
        Reference mask of the unshifted model.
    """

    def setUp(self):
        """Define the filter, DCR parameters, and the bounding box for the tests.
        """
        self.rng = np.random.RandomState(5)
        self.nRandIter = 10  # Number of iterations to repeat each test with random numbers.
        self.dcrNumSubfilters = 3
        self.lambdaEff = 476.31  # Use LSST g band values for the test.
        self.lambdaMin = 405.
        self.lambdaMax = 552.
        self.bufferSize = 5
        xSize = 40
        ySize = 42
        x0 = 12345
        y0 = 67890
        self.bbox = afwGeom.Box2I(afwGeom.Point2I(x0, y0), afwGeom.Extent2I(xSize, ySize))

    def makeTestImages(self, seed=5, nSrc=5, psfSize=2., noiseLevel=5.,
                       detectionSigma=5., sourceSigma=20., fluxRange=2.):
        """Make reproduceable PSF-convolved masked images for testing.

        Parameters
        ----------
        seed : `int`, optional
            Seed value to initialize the random number generator.
        nSrc : `int`, optional
            Number of sources to simulate.
        psfSize : `float`, optional
            Width of the PSF of the simulated sources, in pixels.
        noiseLevel : `float`, optional
            Standard deviation of the noise to add to each pixel.
        detectionSigma : `float`, optional
            Threshold amplitude of the image to set the "DETECTED" mask.
        sourceSigma : `float`, optional
            Average amplitude of the simulated sources,
            relative to ``noiseLevel``
        fluxRange : `float`, optional
            Range in flux amplitude of the simulated sources.

        Returns
        -------
        modelImages : `list` of `lsst.afw.image.Image`
            A list of images, each containing the model for one subfilter
        """
        rng = np.random.RandomState(seed)
        x0, y0 = self.bbox.getBegin()
        xSize, ySize = self.bbox.getDimensions()
        xLoc = rng.rand(nSrc)*(xSize - 2*self.bufferSize) + self.bufferSize + x0
        yLoc = rng.rand(nSrc)*(ySize - 2*self.bufferSize) + self.bufferSize + y0
        modelImages = []

        imageSum = np.zeros((ySize, xSize))
        for subfilter in range(self.dcrNumSubfilters):
            flux = (rng.rand(nSrc)*(fluxRange - 1.) + 1.)*sourceSigma*noiseLevel
            sigmas = [psfSize for src in range(nSrc)]
            coordList = list(zip(xLoc, yLoc, flux, sigmas))
            model = plantSources(self.bbox, 10, 0, coordList, addPoissonNoise=False)
            model.image.array += rng.rand(ySize, xSize)*noiseLevel
            imageSum += model.image.array
            model.mask.addMaskPlane("CLIPPED")
            modelImages.append(model.image)
        maskVals = np.zeros_like(imageSum)
        maskVals[imageSum > detectionSigma*noiseLevel] = afwImage.Mask.getPlaneBitMask('DETECTED')
        model.mask.array[:] = maskVals
        self.mask = model.mask
        return modelImages

    def prepareStats(self):
        """Make a simple statistics object for testing.

        Returns
        -------
        statsCtrl : `lsst.afw.math.StatisticsControl`
            Statistics control object for coaddition.
        """
        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNumSigmaClip(5)
        statsCtrl.setNumIter(3)
        statsCtrl.setNanSafe(True)
        statsCtrl.setWeighted(True)
        statsCtrl.setCalcErrorFromInputVariance(False)
        return statsCtrl

    def makeDummyWcs(self, rotAngle, pixelScale, crval, flipX=True):
        """Make a World Coordinate System object for testing.

        Parameters
        ----------
        rotAngle : `lsst.geom.Angle`
            rotation of the CD matrix, East from North
        pixelScale : `lsst.geom.Angle`
            Pixel scale of the projection.
        crval : `lsst.afw.geom.SpherePoint`
            Coordinates of the reference pixel of the wcs.
        flipX : `bool`, optional
            Flip the direction of increasing Right Ascension.

        Returns
        -------
        `lsst.afw.geom.skyWcs.SkyWcs`
            A wcs that matches the inputs.
        """
        crpix = afwGeom.Box2D(self.bbox).getCenter()
        cdMatrix = afwGeom.makeCdMatrix(scale=pixelScale, orientation=rotAngle, flipX=flipX)
        wcs = afwGeom.makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cdMatrix)
        return wcs

    def makeDummyVisitInfo(self, azimuth, elevation, exposureId=12345, randomizeTime=False):
        """Make a self-consistent visitInfo object for testing.

        Parameters
        ----------
        azimuth : `lsst.geom.Angle`
            Azimuth angle of the simulated observation.
        elevation : `lsst.geom.Angle`
            Elevation angle of the simulated observation.
        exposureId : `int`, optional
            Unique integer identifier for this observation.
        randomizeTime : `bool`, optional
            Add a random offset to the observation time.

        Returns
        -------
        `lsst.afw.image.VisitInfo`
            VisitInfo for the exposure.
        """
        lsstLat = -30.244639*u.degree
        lsstLon = -70.749417*u.degree
        lsstAlt = 2663.*u.m
        lsstTemperature = 20.*u.Celsius
        lsstHumidity = 40.  # in percent
        lsstPressure = 73892.*u.pascal
        loc = EarthLocation(lat=lsstLat,
                            lon=lsstLon,
                            height=lsstAlt)
        airmass = 1.0/np.sin(elevation.asDegrees())

        time = Time(2000., format='decimalyear', scale='tai', location=loc)
        time += 0.5*u.day  # Add half a day since J2000 is defined at noon
        if randomizeTime:
            # Pick a random date and time within a 20-year span
            time += 20*u.year*self.rng.rand()
        altaz = SkyCoord(alt=elevation.asDegrees(), az=azimuth.asDegrees(),
                         unit='deg', obstime=time, frame='altaz', location=loc)
        obsInfo = makeObservationInfo(location=loc,
                                      detector_exposure_id=exposureId,
                                      datetime_begin=time,
                                      datetime_end=time,
                                      boresight_airmass=airmass,
                                      boresight_rotation_angle=Angle(0.*u.degree),
                                      boresight_rotation_coord='sky',
                                      temperature=lsstTemperature,
                                      pressure=lsstPressure,
                                      relative_humidity=lsstHumidity,
                                      tracking_radec=altaz.icrs,
                                      altaz_begin=altaz,
                                      observation_type='science',
                                      )
        makeVisitInfo = MakeRawVisitInfoViaObsInfo()
        visitInfo = makeVisitInfo.observationInfo2visitInfo(obsInfo=obsInfo)
        return visitInfo

    def testDummyVisitInfo(self):
        """Verify the implementation of the visitInfo used for tests.
        """
        azimuth = 0*degrees
        for testIter in range(self.nRandIter):
            # Restrict to 45 < elevation < 85 degrees
            elevation = (45. + self.rng.rand()*40.)*degrees
            visitInfo = self.makeDummyVisitInfo(azimuth, elevation)
            dec = visitInfo.getBoresightRaDec().getLatitude()
            lat = visitInfo.getObservatory().getLatitude()
            # An observation made with azimuth=0 should be pointed to the North
            # So the RA should be North of the telescope's latitude
            self.assertGreater(dec.asDegrees(), lat.asDegrees())

            # The hour angle should be zero for azimuth=0
            HA = visitInfo.getBoresightHourAngle()
            refHA = 0.*degrees
            self.assertAnglesAlmostEqual(HA, refHA, maxDiff=coordinateTolerance)
            # If the observation is North of the telescope's latitude, the
            # direction to zenith should be along the -y axis
            # with a parallactic angle of 180 degrees
            parAngle = visitInfo.getBoresightParAngle()
            refParAngle = 180.*degrees
            self.assertAnglesAlmostEqual(parAngle, refParAngle, maxDiff=coordinateTolerance)

    def testDcrCalculation(self):
        """Test that the shift in pixels due to DCR is consistently computed.

        The shift is compared to pre-computed values.
        """
        dcrNumSubfilters = 3
        afwImageUtils.defineFilter("gTest", self.lambdaEff,
                                   lambdaMin=self.lambdaMin, lambdaMax=self.lambdaMax)
        filterInfo = afwImage.Filter("gTest")
        rotAngle = 0.*degrees
        azimuth = 30.*degrees
        elevation = 65.*degrees
        pixelScale = 0.2*arcseconds
        visitInfo = self.makeDummyVisitInfo(azimuth, elevation)
        wcs = self.makeDummyWcs(rotAngle, pixelScale, crval=visitInfo.getBoresightRaDec())
        dcrShift = calculateDcr(visitInfo, wcs, filterInfo, dcrNumSubfilters)
        # Compare to precomputed values.
        refShift = [(-0.5575567724366292, -0.2704095599533037),
                    (0.001961910992342903, 0.000951507567181944),
                    (0.40402552599550073, 0.19594841296051665)]
        for shiftOld, shiftNew in zip(refShift, dcrShift):
            self.assertFloatsAlmostEqual(shiftOld[1], shiftNew[1], rtol=1e-6, atol=1e-8)
            self.assertFloatsAlmostEqual(shiftOld[0], shiftNew[0], rtol=1e-6, atol=1e-8)

    def testCoordinateTransformDcrCalculation(self):
        """Check the DCR calculation using astropy coordinate transformations.

        Astmospheric refraction causes sources to appear closer to zenith than
        they really are. An alternate calculation of the shift due to DCR is to
        transform the pixel coordinates to altitude and azimuth, add the DCR
        amplitude to the altitude, and transform back to pixel coordinates.
        """
        dcrNumSubfilters = 3
        afwImageUtils.defineFilter("gTest", self.lambdaEff,
                                   lambdaMin=self.lambdaMin, lambdaMax=self.lambdaMax)
        filterInfo = afwImage.Filter("gTest")
        lambdaEff = filterInfo.getFilterProperty().getLambdaEff()
        pixelScale = 0.2*arcseconds

        for testIter in range(self.nRandIter):
            rotAngle = 360.*self.rng.rand()*degrees
            azimuth = 360.*self.rng.rand()*degrees
            elevation = (45. + self.rng.rand()*40.)*degrees  # Restrict to 45 < elevation < 85 degrees
            visitInfo = self.makeDummyVisitInfo(azimuth, elevation)
            loc = EarthLocation(lat=visitInfo.getObservatory().getLatitude().asDegrees()*u.degree,
                                lon=visitInfo.getObservatory().getLongitude().asDegrees()*u.degree,
                                height=visitInfo.getObservatory().getElevation()*u.m)
            wcs = self.makeDummyWcs(rotAngle, pixelScale, crval=visitInfo.getBoresightRaDec())
            date = visitInfo.getDate()
            time = Time(date.get(date.MJD, date.TAI), format='mjd', location=loc, scale='tai')
            altaz = SkyCoord(alt=elevation.asDegrees(), az=azimuth.asDegrees(),
                             unit='deg', obstime=time, frame='altaz', location=loc)
            # The DCR calculations are performed at the boresight
            ra0 = altaz.icrs.ra.degree*degrees
            dec0 = altaz.icrs.dec.degree*degrees
            x0, y0 = wcs.skyToPixel(afwGeom.SpherePoint(ra0, dec0))
            dcrShifts = calculateDcr(visitInfo, wcs, filterInfo, dcrNumSubfilters)
            refShifts = []
            # We divide the filter into "subfilters" with with the full wavelength range
            # divided into equal sub-ranges.
            for wl0, wl1 in wavelengthGenerator(filterInfo, dcrNumSubfilters):
                # Note that diffRefractAmp can be negative,
                # since it is relative to the midpoint of the full band
                diffRefractAmp0 = differentialRefraction(wavelength=wl0, wavelengthRef=lambdaEff,
                                                         elevation=elevation,
                                                         observatory=visitInfo.getObservatory(),
                                                         weather=visitInfo.getWeather())
                diffRefractAmp1 = differentialRefraction(wavelength=wl1, wavelengthRef=lambdaEff,
                                                         elevation=elevation,
                                                         observatory=visitInfo.getObservatory(),
                                                         weather=visitInfo.getWeather())
                diffRefractAmp = (diffRefractAmp0 + diffRefractAmp1)/2.

                elevation1 = elevation + diffRefractAmp
                altaz = SkyCoord(alt=elevation1.asDegrees(), az=azimuth.asDegrees(),
                                 unit='deg', obstime=time, frame='altaz', location=loc)
                ra1 = altaz.icrs.ra.degree*degrees
                dec1 = altaz.icrs.dec.degree*degrees
                x1, y1 = wcs.skyToPixel(afwGeom.SpherePoint(ra1, dec1))
                refShifts.append((y1-y0, x1-x0))
            for refShift, dcrShift in zip(refShifts, dcrShifts):
                # Use a fairly loose tolerance, since 1% of a pixel is good enough agreement.
                self.assertFloatsAlmostEqual(refShift[1], dcrShift[1], rtol=1e-2, atol=1e-2)
                self.assertFloatsAlmostEqual(refShift[0], dcrShift[0], rtol=1e-2, atol=1e-2)

    def testDcrSubfilterOrder(self):
        """Test that the bluest subfilter always has the largest DCR amplitude.
        """
        dcrNumSubfilters = 3
        afwImageUtils.defineFilter("gTest", self.lambdaEff,
                                   lambdaMin=self.lambdaMin, lambdaMax=self.lambdaMax)
        filterInfo = afwImage.Filter("gTest")
        pixelScale = 0.2*arcseconds
        for testIter in range(self.nRandIter):
            rotAngle = 360.*self.rng.rand()*degrees
            azimuth = 360.*self.rng.rand()*degrees
            elevation = (45. + self.rng.rand()*40.)*degrees  # Restrict to 45 < elevation < 85 degrees
            visitInfo = self.makeDummyVisitInfo(azimuth, elevation)
            wcs = self.makeDummyWcs(rotAngle, pixelScale, crval=visitInfo.getBoresightRaDec())
            dcrShift = calculateDcr(visitInfo, wcs, filterInfo, dcrNumSubfilters)
            # First check that the blue subfilter amplitude is greater than the red subfilter
            rotation = calculateImageParallacticAngle(visitInfo, wcs).asRadians()
            ampShift = [dcr[1]*np.sin(rotation) + dcr[0]*np.cos(rotation) for dcr in dcrShift]
            self.assertGreater(ampShift[0], 0.)  # The blue subfilter should be shifted towards zenith
            self.assertLess(ampShift[2], 0.)  # The red subfilter should be shifted away from zenith
            # The absolute amplitude of the blue subfilter should also
            # be greater than that of the red subfilter
            self.assertGreater(np.abs(ampShift[0]), np.abs(ampShift[2]))

    def testApplyDcr(self):
        """Test that the image transformation reduces to a simple shift.
        """
        dxVals = [-2, 1, 0, 1, 2]
        dyVals = [-2, 1, 0, 1, 2]
        x0 = 13
        y0 = 27
        inputImage = afwImage.MaskedImageF(self.bbox)
        image = inputImage.image.array
        image[y0, x0] = 1.
        for dx in dxVals:
            for dy in dyVals:
                shift = (dy, dx)
                shiftedImage = applyDcr(image, shift, useInverse=False)
                # Create a blank reference image, and add the fake point source at the shifted location.
                refImage = afwImage.MaskedImageF(self.bbox)
                refImage.image.array[y0 + dy, x0 + dx] = 1.
                self.assertFloatsAlmostEqual(shiftedImage, refImage.image.array, rtol=1e-12, atol=1e-12)

    def testRotationAngle(self):
        """Test that the sky rotation angle is consistently computed.

        The rotation is compared to pre-computed values.
        """
        cdRotAngle = 0.*degrees
        azimuth = 130.*degrees
        elevation = 70.*degrees
        pixelScale = 0.2*arcseconds
        visitInfo = self.makeDummyVisitInfo(azimuth, elevation)
        wcs = self.makeDummyWcs(cdRotAngle, pixelScale, crval=visitInfo.getBoresightRaDec())
        rotAngle = calculateImageParallacticAngle(visitInfo, wcs)
        refAngle = -1.0848040464064805*radians
        self.assertAnglesAlmostEqual(refAngle, rotAngle)

    def testRotationSouthZero(self):
        """Test that an observation pointed due South has zero rotation angle.

        An observation pointed South and on the meridian should have zenith
        directly to the North, and a parallactic angle of zero.
        """
        refAngle = 0.*degrees
        azimuth = 180.*degrees  # Telescope is pointed South
        pixelScale = 0.2*arcseconds
        for testIter in range(self.nRandIter):
            # Any additional arbitrary rotation should fall out of the calculation
            cdRotAngle = 360*self.rng.rand()*degrees
            elevation = (45. + self.rng.rand()*40.)*degrees  # Restrict to 45 < elevation < 85 degrees
            visitInfo = self.makeDummyVisitInfo(azimuth, elevation)
            wcs = self.makeDummyWcs(cdRotAngle, pixelScale, crval=visitInfo.getBoresightRaDec(), flipX=True)
            rotAngle = calculateImageParallacticAngle(visitInfo, wcs)
            self.assertAnglesAlmostEqual(refAngle - cdRotAngle, rotAngle, maxDiff=coordinateTolerance)

    def testRotationFlipped(self):
        """Check the interpretation of rotations in the WCS.
        """
        doFlip = [False, True]
        for testIter in range(self.nRandIter):
            # Any additional arbitrary rotation should fall out of the calculation
            cdRotAngle = 360*self.rng.rand()*degrees
            # Make the telescope be pointed South, so that the parallactic angle is zero.
            azimuth = 180.*degrees
            elevation = (45. + self.rng.rand()*40.)*degrees  # Restrict to 45 < elevation < 85 degrees
            pixelScale = 0.2*arcseconds
            visitInfo = self.makeDummyVisitInfo(azimuth, elevation)
            for flip in doFlip:
                wcs = self.makeDummyWcs(cdRotAngle, pixelScale,
                                        crval=visitInfo.getBoresightRaDec(),
                                        flipX=flip)
                rotAngle = calculateImageParallacticAngle(visitInfo, wcs)
                if flip:
                    rotAngle *= -1
                self.assertAnglesAlmostEqual(cdRotAngle, rotAngle, maxDiff=coordinateTolerance)

    def testConditionDcrModelNoChange(self):
        """Conditioning should not change the model if it equals the reference.
        """
        modelImages = self.makeTestImages()
        dcrModels = DcrModel(modelImages=modelImages, mask=self.mask)
        newModels = [model.clone() for model in dcrModels]
        dcrModels.conditionDcrModel(newModels, self.bbox, gain=1.)
        for refModel, newModel in zip(dcrModels, newModels):
            self.assertFloatsAlmostEqual(refModel.array, newModel.array)

    def testConditionDcrModelNoChangeHighGain(self):
        """Conditioning should not change the model if it equals the reference.
        """
        modelImages = self.makeTestImages()
        dcrModels = DcrModel(modelImages=modelImages, mask=self.mask)
        newModels = [model.clone() for model in dcrModels]
        dcrModels.conditionDcrModel(newModels, self.bbox, gain=3.)
        for refModel, newModel in zip(dcrModels, newModels):
            self.assertFloatsAlmostEqual(refModel.array, newModel.array)

    def testConditionDcrModelWithChange(self):
        """Verify conditioning when the model changes by a known amount.
        """
        modelImages = self.makeTestImages()
        dcrModels = DcrModel(modelImages=modelImages, mask=self.mask)
        newModels = [model.clone() for model in dcrModels]
        for model in newModels:
            model.array[:] *= 3.
        dcrModels.conditionDcrModel(newModels, self.bbox, gain=1.)
        for refModel, newModel in zip(dcrModels, newModels):
            refModel.array[:] *= 2.
            self.assertFloatsAlmostEqual(refModel.array, newModel.array)

    def testRegularizationSmallClamp(self):
        """Test that large variations between model planes are reduced.

        This also tests that noise-like pixels are not regularized.
        """
        clampFrequency = 2
        regularizationWidth = 2
        fluxRange = 10.
        modelImages = self.makeTestImages(fluxRange=fluxRange)
        dcrModels = DcrModel(modelImages=modelImages, mask=self.mask)
        newModels = [model.clone() for model in dcrModels]
        templateImage = dcrModels.getReferenceImage(self.bbox)

        statsCtrl = self.prepareStats()
        dcrModels.regularizeModelFreq(newModels, self.bbox, statsCtrl, clampFrequency, regularizationWidth)
        for model, refModel in zip(newModels, dcrModels):
            # Make sure the test parameters do reduce the outliers
            self.assertGreater(np.max(refModel.array - templateImage),
                               np.max(model.array - templateImage))
            highThreshold = templateImage*clampFrequency
            highPix = model.array > highThreshold
            highPix = ndimage.morphology.binary_opening(highPix, iterations=regularizationWidth)
            self.assertFalse(np.all(highPix))
            lowThreshold = templateImage/clampFrequency
            lowPix = model.array < lowThreshold
            lowPix = ndimage.morphology.binary_opening(lowPix, iterations=regularizationWidth)
            self.assertFalse(np.all(lowPix))

    def testRegularizationSidelobes(self):
        """Test that artificial chromatic sidelobes are suppressed.
        """
        clampFrequency = 2.
        regularizationWidth = 2
        noiseLevel = 0.1
        sourceAmplitude = 100.
        modelImages = self.makeTestImages(seed=5, nSrc=5, psfSize=3., noiseLevel=noiseLevel,
                                          detectionSigma=5., sourceSigma=sourceAmplitude, fluxRange=2.)
        templateImage = np.mean([model.array for model in modelImages], axis=0)
        sidelobeImages = self.makeTestImages(seed=5, nSrc=5, psfSize=1.5, noiseLevel=noiseLevel/10.,
                                             detectionSigma=5., sourceSigma=sourceAmplitude*5., fluxRange=2.)
        statsCtrl = self.prepareStats()
        signList = [-1., 0., 1.]
        sidelobeShift = (0., 4.)
        for model, sidelobe, sign in zip(modelImages, sidelobeImages, signList):
            sidelobe.array *= sign
            model.array += applyDcr(sidelobe.array, sidelobeShift, useInverse=False)
            model.array += applyDcr(sidelobe.array, sidelobeShift, useInverse=True)

        dcrModels = DcrModel(modelImages=modelImages, mask=self.mask)
        refModels = [dcrModels[subfilter].clone() for subfilter in range(self.dcrNumSubfilters)]

        dcrModels.regularizeModelFreq(modelImages, self.bbox, statsCtrl, clampFrequency,
                                      regularizationWidth=regularizationWidth)
        for model, refModel, sign in zip(modelImages, refModels, signList):
            # Make sure the test parameters do reduce the outliers
            self.assertGreater(np.sum(np.abs(refModel.array - templateImage)),
                               np.sum(np.abs(model.array - templateImage)))

    def testRegularizeModelIter(self):
        """Test that large amplitude changes between iterations are restricted.

        This also tests that noise-like pixels are not regularized.
        """
        modelClampFactor = 2.
        regularizationWidth = 2
        subfilter = 0
        dcrModels = DcrModel(modelImages=self.makeTestImages())
        oldModel = dcrModels[0]
        xSize, ySize = self.bbox.getDimensions()
        newModel = oldModel.clone()
        newModel.array[:] += self.rng.rand(ySize, xSize)*np.max(oldModel.array)
        newModelRef = newModel.clone()

        dcrModels.regularizeModelIter(subfilter, newModel, self.bbox, modelClampFactor, regularizationWidth)

        # Make sure the test parameters do reduce the outliers
        self.assertGreater(np.max(newModelRef.array),
                           np.max(newModel.array - oldModel.array))
        # Check that all of the outliers are clipped
        highThreshold = oldModel.array*modelClampFactor
        highPix = newModel.array > highThreshold
        highPix = ndimage.morphology.binary_opening(highPix, iterations=regularizationWidth)
        self.assertFalse(np.all(highPix))
        lowThreshold = oldModel.array/modelClampFactor
        lowPix = newModel.array < lowThreshold
        lowPix = ndimage.morphology.binary_opening(lowPix, iterations=regularizationWidth)
        self.assertFalse(np.all(lowPix))

    def testIterateModel(self):
        """Test that the DcrModel is iterable, and has the right values.
        """
        testModels = self.makeTestImages()
        refVals = [np.sum(model.array) for model in testModels]
        dcrModels = DcrModel(modelImages=testModels)
        for refVal, model in zip(refVals, dcrModels):
            self.assertFloatsEqual(refVal, np.sum(model.array))
        # Negative indices are allowed, so check that those return models from the end.
        self.assertFloatsEqual(refVals[-1], np.sum(dcrModels[-1].array))


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
