# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import os
import unittest
from unittest import mock
import requests

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.geom
from lsst.ip.diffim import detectAndMeasure, subtractImages
from lsst.afw.table import IdFactory
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.meas.algorithms as measAlg
from lsst.pipe.base import InvalidQuantumError, UpstreamFailureNoWorkFound, AlgorithmError
import lsst.utils.tests
import lsst.meas.base.tests
import lsst.daf.base as dafBase
from lsst.afw.coord import Observatory, Weather
import lsst.geom as geom
import lsst.pex.config as pexConfig

from utils import makeTestImage, checkMask


class DetectAndMeasureTestBase:

    def _check_diaSource(self, refSources, diaSource, refIds=None,
                         matchDistance=1., scale=1., usePsfFlux=True,
                         rtol=0.025, atol=None):
        """Match a diaSource with a source in a reference catalog
        and compare properties.

        Parameters
        ----------
        refSources : `lsst.afw.table.SourceCatalog`
            The reference catalog.
        diaSource : `lsst.afw.table.SourceRecord`
            The new diaSource to match to the reference catalog.
        refIds : `list` of `int`, optional
            Source IDs of previously associated diaSources.
        matchDistance : `float`, optional
            Maximum distance allowed between the detected and reference source
            locations, in pixels.
        scale : `float`, optional
            Optional factor to scale the flux by before performing the test.
        usePsfFlux : `bool`, optional
            If set, test the PsfInstFlux field, otherwise use ApInstFlux.
        rtol : `float`, optional
            Relative tolerance of the flux value test.
        atol : `float`, optional
            Absolute tolerance of the flux value test.
        """
        if diaSource["sky_source"]:
            return
        distance = np.sqrt((diaSource.getX() - refSources.getX())**2
                           + (diaSource.getY() - refSources.getY())**2)
        self.assertLess(min(distance), matchDistance)
        src = refSources[np.argmin(distance)]
        if refIds is not None:
            # Check that the same source was not previously associated
            self.assertNotIn(src.getId(), refIds)
            refIds.append(src.getId())
        if atol is None:
            atol = rtol*src.getPsfInstFlux() if usePsfFlux else rtol*src.getApInstFlux()
        if usePsfFlux:
            self.assertFloatsAlmostEqual(src.getPsfInstFlux()*scale, diaSource.getPsfInstFlux(),
                                         rtol=rtol, atol=atol)
        else:
            self.assertFloatsAlmostEqual(src.getApInstFlux()*scale, diaSource.getApInstFlux(),
                                         rtol=rtol, atol=atol)

    def _check_values(self, values, minValue=None, maxValue=None):
        """Verify that an array has finite values, and optionally that they are
        within specified minimum and maximum bounds.

        Parameters
        ----------
        values : `numpy.ndarray`
            Array of values to check.
        minValue : `float`, optional
            Minimum allowable value.
        maxValue : `float`, optional
            Maximum allowable value.
        """
        self.assertTrue(np.all(np.isfinite(values)))
        if minValue is not None:
            self.assertTrue(np.all(values >= minValue))
        if maxValue is not None:
            self.assertTrue(np.all(values <= maxValue))

    def _setup_detection(self, doSkySources=True, nSkySources=5,
                         doSubtractBackground=False, run_sattle=False, **kwargs):
        """Setup and configure the detection and measurement PipelineTask.

        Parameters
        ----------
        doSkySources : `bool`, optional
            Generate sky sources.
        nSkySources : `int`, optional
            The number of sky sources to add in isolated background regions.
        **kwargs
            Any additional config parameters to set.

        Returns
        -------
        `lsst.pipe.base.PipelineTask`
            The configured Task to use for detection and measurement.
        """
        config = self.detectionTask.ConfigClass()
        config.doSkySources = doSkySources
        if doSkySources:
            config.skySources.nSources = nSkySources
        config.update(**kwargs)

        config.run_sattle = run_sattle

        # Make a realistic id generator so that output catalog ids are useful.
        dataId = lsst.daf.butler.DataCoordinate.standardize(
            instrument="I",
            visit=42,
            detector=12,
            universe=lsst.daf.butler.DimensionUniverse(),
        )
        config.idGenerator.packer.name = "observation"
        config.idGenerator.packer["observation"].n_observations = 10000
        config.idGenerator.packer["observation"].n_detectors = 99
        config.idGenerator.n_releases = 8
        config.idGenerator.release_id = 2
        config.doSubtractBackground = doSubtractBackground
        self.idGenerator = config.idGenerator.apply(dataId)

        return self.detectionTask(config=config)


class DetectAndMeasureTest(DetectAndMeasureTestBase, lsst.utils.tests.TestCase):
    detectionTask = detectAndMeasure.DetectAndMeasureTask

    def test_detection_xy0(self):
        """Basic functionality test with non-zero x0 and y0.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "x0": 12345, "y0": 67890}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo())
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()

        # Configure the detection Task
        # Set the subtraction residual metric threshold high, since the template
        # is not actually subtracted from the science image.
        detectionTask = self._setup_detection(doDeblend=True, badSubtractionRatioThreshold=1.,
                                              doSkySources=False)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference, sources,
                                   idFactory=self.idGenerator.make_table_id_factory())
        subtractedMeasuredExposure = output.subtractedMeasuredExposure

        # Catalog ids should be very large from this id generator.
        self.assertTrue(all(output.diaSources['id'] > 1000000000))
        self.assertImagesEqual(subtractedMeasuredExposure.image, difference.image)

        # all of the sources should have been detected
        self.assertEqual(len(output.diaSources), len(sources))
        # no sources should be flagged as negative
        self.assertEqual(len(~output.diaSources["is_negative"]), len(output.diaSources))
        refIds = []
        for source in sources:
            self._check_diaSource(output.diaSources, source, refIds=refIds)

    def test_raise_bad_psf(self):
        """Detection should raise if the PSF width is NaN
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "x0": 12345, "y0": 67890}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()

        psfShape = difference.getPsf().computeBBox(lsst.geom.Point2D(0, 0)).getDimensions()
        psfcI = afwImage.ImageD(lsst.geom.Extent2I(psfShape[1], psfShape[0]))
        psfNew = np.zeros(psfShape)
        psfNew[0, :] = 1
        psfcI.array = psfNew
        psfcK = afwMath.FixedKernel(psfcI)
        psf = measAlg.KernelPsf(psfcK)
        difference.setPsf(psf)

        # Configure the detection Task
        detectionTask = self._setup_detection()

        # Run detection and check the results
        with self.assertRaises(UpstreamFailureNoWorkFound):
            detectionTask.run(science, matchedTemplate, difference, sources)

    def test_measurements_finite(self):
        """Measured fluxes and centroids should always be finite.
        """
        columnNames = ["coord_ra", "coord_dec", "ip_diffim_forced_PsfFlux_instFlux"]

        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        xSize = 256
        ySize = 256
        kwargs = {"psfSize": 2.4, "x0": 0, "y0": 0,
                  "xSize": xSize, "ySize": ySize}
        science, sources = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel, noiseSeed=6,
                                         nSrc=1, **kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo())
        matchedTemplate, _ = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel/4, noiseSeed=7,
                                           nSrc=1, **kwargs)
        rng = np.random.RandomState(3)
        xLoc = np.arange(-5, xSize+5, 10)
        rng.shuffle(xLoc)
        yLoc = np.arange(-5, ySize+5, 10)
        rng.shuffle(yLoc)
        transients, transientSources = makeTestImage(seed=transientSeed,
                                                     nSrc=len(xLoc), fluxLevel=1000.,
                                                     noiseLevel=noiseLevel, noiseSeed=8,
                                                     xLoc=xLoc, yLoc=yLoc,
                                                     **kwargs)
        difference = science.clone()
        difference.maskedImage -= matchedTemplate.maskedImage
        difference.maskedImage += transients.maskedImage

        # Configure the detection Task
        detectionTask = self._setup_detection(doForcedMeasurement=True)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference, sources)

        for column in columnNames:
            self._check_values(output.diaSources[column])
        self._check_values(output.diaSources.getX(), minValue=0, maxValue=xSize)
        self._check_values(output.diaSources.getY(), minValue=0, maxValue=ySize)
        self._check_values(output.diaSources.getPsfInstFlux())

    def test_raise_config_schema_mismatch(self):
        """Check that sources with specified flags are removed from the catalog.
        """
        # Configure the detection Task, and and set a config that is not in the schema
        with self.assertRaises(InvalidQuantumError):
            self._setup_detection(badSourceFlags=["Bogus_flag_42"])

    def test_remove_unphysical(self):
        """Check that sources with specified flags are removed from the catalog.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        xSize = 256
        ySize = 256
        kwargs = {"psfSize": 2.4, "xSize": xSize, "ySize": ySize}
        science, sources = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel, noiseSeed=6,
                                         nSrc=1, **kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo())
        matchedTemplate, _ = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel/4, noiseSeed=7,
                                           nSrc=1, **kwargs)
        transients, transientSources = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=8, nSrc=1, **kwargs)
        science.maskedImage += transients.maskedImage
        difference = science.clone()
        bbox = difference.getBBox()
        difference.maskedImage -= matchedTemplate.maskedImage

        # Configure the detection Task, and remove unphysical sources
        detectionTask = self._setup_detection(doForcedMeasurement=False, doSkySources=True, nSkySources=20,
                                              badSourceFlags=["base_PixelFlags_flag_offimage", ])

        # Run detection and check the results
        diaSources = detectionTask.run(science, matchedTemplate, difference, sources).diaSources
        badDiaSrcDoRemove = ~bbox.contains(diaSources.getX(), diaSources.getY())
        nBadDoRemove = np.count_nonzero(badDiaSrcDoRemove)
        # Verify that all sources are physical
        self.assertEqual(nBadDoRemove, 0)
        # Set a few centroids outside the image bounding box
        nSetBad = 5
        for src in diaSources[0: nSetBad]:
            src["slot_Centroid_x"] += xSize
            src["slot_Centroid_y"] += ySize
            src["base_PixelFlags_flag_offimage"] = True
        # Verify that these sources are outside the image
        badDiaSrc = ~bbox.contains(diaSources.getX(), diaSources.getY())
        nBad = np.count_nonzero(badDiaSrc)
        self.assertEqual(nBad, nSetBad)
        diaSourcesNoBad = detectionTask._removeBadSources(diaSources)
        badDiaSrcNoBad = ~bbox.contains(diaSourcesNoBad.getX(), diaSourcesNoBad.getY())

        # Verify that no sources outside the image bounding box remain
        self.assertEqual(np.count_nonzero(badDiaSrcNoBad), 0)
        self.assertEqual(len(diaSourcesNoBad), len(diaSources) - nSetBad)

    def test_detect_transients(self):
        """Run detection on a difference image containing transients.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        scienceBase, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)

        # Configure the detection Task
        detectionTask = self._setup_detection(doMerge=False, doSkySources=True)
        kwargs["seed"] = transientSeed
        kwargs["nSrc"] = 10
        kwargs["fluxLevel"] = 1000

        # Run detection and check the results
        def _detection_wrapper(positive=True):
            transients, transientSources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=8, **kwargs)
            science = scienceBase.clone()
            if positive:
                science.maskedImage += transients.maskedImage
            else:
                science.maskedImage -= transients.maskedImage
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            # NOTE: NoiseReplacer (run by forcedMeasurement) can modify the
            # science image if we've e.g. removed parents post-deblending.
            # Pass a clone of the science image, so that it doesn't disrupt
            # later tests.
            output = detectionTask.run(science, matchedTemplate, difference, sources)
            refIds = []
            scale = 1. if positive else -1.
            for diaSource in output.diaSources:
                self._check_diaSource(transientSources, diaSource, refIds=refIds, scale=scale)
        _detection_wrapper(positive=True)
        _detection_wrapper(positive=False)

    def test_detect_transients_with_background(self):
        """Run detection on a difference image containing transients and a background.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        fluxLevel = 500
        xSize = 512
        ySize = 512
        x0 = 123
        y0 = 456
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel,
                  "xSize": xSize, "ySize": ySize, "x0": x0, "y0": y0}
        params = [2.2, 2.1, 2.0, 1.2, 1.1, 1.0]

        bbox2D = lsst.geom.Box2D(lsst.geom.Point2D(x0, y0), lsst.geom.Extent2D(xSize, ySize))
        background_model = afwMath.Chebyshev1Function2D(params, bbox2D)
        scienceBase, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6,
                                             background=background_model, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)

        # Configure the detection Task
        detectionTask = self._setup_detection(doMerge=False, doSubtractBackground=True)
        kwargs["seed"] = transientSeed
        kwargs["nSrc"] = 10
        kwargs["fluxLevel"] = 1000

        # Run detection and check the results
        def _detection_wrapper(positive=True):
            transients, transientSources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=8, **kwargs)
            science = scienceBase.clone()
            if positive:
                science.maskedImage += transients.maskedImage
            else:
                science.maskedImage -= transients.maskedImage
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            # NOTE: NoiseReplacer (run by forcedMeasurement) can modify the
            # science image if we've e.g. removed parents post-deblending.
            # Pass a clone of the science image, so that it doesn't disrupt
            # later tests.
            output = detectionTask.run(science, matchedTemplate, difference, sources)
            refIds = []
            scale = 1. if positive else -1.
            for transient in transientSources:
                self._check_diaSource(output.diaSources, transient, refIds=refIds, scale=scale)
        _detection_wrapper(positive=True)
        _detection_wrapper(positive=False)

    def test_missing_mask_planes(self):
        """Check that detection runs with missing mask planes.
        """
        # Set up the simulated images
        noiseLevel = 1.
        fluxLevel = 500
        kwargs = {"psfSize": 2.4, "fluxLevel": fluxLevel, "addMaskPlanes": []}
        # Use different seeds for the science and template so every source is a diaSource
        science, sources = makeTestImage(seed=5, noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo())
        matchedTemplate, _ = makeTestImage(seed=6, noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)

        difference = science.clone()
        difference.maskedImage -= matchedTemplate.maskedImage
        detectionTask = self._setup_detection(raiseOnBadSubtractionRatio=False, raiseOnNoDiaSources=False)

        # Verify that detection runs without errors
        detectionTask.run(science, matchedTemplate, difference, sources)

    def test_detect_dipoles(self):
        """Run detection on a difference image containing dipoles.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 1000
        fluxRange = 1.5
        nSources = 10
        offset = 1
        xSize = 300
        ySize = 300
        kernelSize = 32
        # Avoid placing sources near the edge for this test, so that we can
        # easily check that the correct number of sources are detected.
        templateBorderSize = kernelSize//2
        dipoleFlag = "ip_diffim_DipoleFit_classification"
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "fluxRange": fluxRange,
                  "nSrc": nSources, "templateBorderSize": templateBorderSize, "kernelSize": kernelSize,
                  "xSize": xSize, "ySize": ySize}
        dipoleFlag = "ip_diffim_DipoleFit_classification"
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo())
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()
        matchedTemplate.image.array[...] = np.roll(matchedTemplate.image.array[...], offset, axis=0)
        matchedTemplate.variance.array[...] = np.roll(matchedTemplate.variance.array[...], offset, axis=0)
        matchedTemplate.mask.array[...] = np.roll(matchedTemplate.mask.array[...], offset, axis=0)
        difference.maskedImage -= matchedTemplate.maskedImage[science.getBBox()]

        # Need a higher residual metric threshold, since we are purposely
        # creating poor subtractions.
        detectionTask = self._setup_detection(doMerge=True, badSubtractionRatioThreshold=0.3)
        output = detectionTask.run(science, matchedTemplate, difference, sources)
        diaSources = output.diaSources[~output.diaSources["sky_source"]].copy(deep=True)
        self.assertEqual(len(diaSources), len(sources))
        # no sources should be flagged as negative
        self.assertEqual(len(~diaSources["is_negative"]), len(diaSources))
        refIds = []
        for diaSource in diaSources:
            if diaSource[dipoleFlag]:
                self._check_diaSource(sources, diaSource, refIds=refIds, scale=0,
                                      rtol=0.05, atol=None, usePsfFlux=False)
                self.assertFloatsAlmostEqual(diaSource["ip_diffim_DipoleFit_orientation"],
                                             -np.pi / 2, atol=2.)
                self.assertFloatsAlmostEqual(diaSource["ip_diffim_DipoleFit_separation"], offset, rtol=0.1)
            else:
                raise ValueError("DiaSource with ID %s is not a dipole!", diaSource.getId())

    def test_sky_sources(self):
        """Add sky sources and check that they are sufficiently far from other
        sources and have negligible flux.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        transientFluxLevel = 1000.
        transientFluxRange = 1.5
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo())
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        transients, transientSources = makeTestImage(seed=transientSeed, psfSize=2.4,
                                                     nSrc=10, fluxLevel=transientFluxLevel,
                                                     fluxRange=transientFluxRange,
                                                     noiseLevel=noiseLevel, noiseSeed=8)
        difference = science.clone()
        difference.maskedImage -= matchedTemplate.maskedImage
        difference.maskedImage += transients.maskedImage
        kernelWidth = np.max(science.psf.getKernel().getDimensions())//2

        # Configure the detection Task
        detectionTask = self._setup_detection(doSkySources=True)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference, sources,
                                   idFactory=self.idGenerator.make_table_id_factory())
        skySources = output.diaSources[output.diaSources["sky_source"]]
        self.assertEqual(len(skySources), detectionTask.config.skySources.nSources)
        for skySource in skySources:
            # Disable the sky_source flag to enable checks.
            skySource["sky_source"] = False
            # The sky sources should not be close to any other source
            with self.assertRaises(AssertionError):
                self._check_diaSource(transientSources, skySource, matchDistance=kernelWidth)
            with self.assertRaises(AssertionError):
                self._check_diaSource(sources, skySource, matchDistance=kernelWidth)
            # The sky sources should have low flux levels.
            self._check_diaSource(transientSources, skySource, matchDistance=1000, scale=0.,
                                  atol=np.sqrt(transientFluxRange*transientFluxLevel))

        # Catalog ids should be very large from this id generator.
        self.assertTrue(all(output.diaSources['id'] > 1000000000))

    def test_exclude_mask_detections(self):
        """Sources with certain bad mask planes set should not be detected.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        fluxLevel = 500
        radius = 2
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo())
        detector = DetectorWrapper(numAmps=1).detector
        science.setDetector(detector)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)

        # Configure the detection Task
        detectionTask = self._setup_detection()
        excludeMaskPlanes = detectionTask.config.detection.excludeMaskPlanes
        nBad = len(excludeMaskPlanes)
        self.assertGreater(nBad, 0)
        kwargs["seed"] = transientSeed
        kwargs["nSrc"] = nBad
        kwargs["fluxLevel"] = 1000

        # Run detection and check the results
        def _detection_wrapper(setFlags=True):
            transients, transientSources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=8, **kwargs)
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            difference.maskedImage += transients.maskedImage
            if setFlags:
                for src, badMask in zip(transientSources, excludeMaskPlanes):
                    srcX = int(src.getX())
                    srcY = int(src.getY())
                    srcBbox = lsst.geom.Box2I(lsst.geom.Point2I(srcX - radius, srcY - radius),
                                              lsst.geom.Extent2I(2*radius + 1, 2*radius + 1))
                    difference[srcBbox].mask.array |= lsst.afw.image.Mask.getPlaneBitMask(badMask)
            if setFlags:
                with self.assertRaises(AlgorithmError):
                    output = detectionTask.run(science, matchedTemplate, difference, sources)
                return
            else:
                output = detectionTask.run(science, matchedTemplate, difference, sources)
                refIds = []
                goodSrcFlags = checkMask(difference.mask, transientSources, excludeMaskPlanes)
                self.assertEqual(np.sum(~goodSrcFlags), 0)
                for diaSource, goodSrcFlag in zip(output.diaSources, goodSrcFlags):
                    if ~goodSrcFlag:
                        with self.assertRaises(AssertionError):
                            self._check_diaSource(transientSources, diaSource, refIds=refIds)
                    else:
                        self._check_diaSource(transientSources, diaSource, refIds=refIds)
        _detection_wrapper(setFlags=False)
        _detection_wrapper(setFlags=True)

    def test_fake_mask_plane_propagation(self):
        """Test that we have the mask planes related to fakes in diffim images.
        This is testing method called updateMasks
        """
        xSize = 256
        ySize = 256
        science, sources = makeTestImage(psfSize=2.4, xSize=xSize, ySize=ySize, doApplyCalibration=True)
        science_fake_img, science_fake_sources = makeTestImage(
            psfSize=2.4, xSize=xSize, ySize=ySize, seed=5, nSrc=3, noiseLevel=0.25, fluxRange=1
        )
        template, _ = makeTestImage(psfSize=2.4, xSize=xSize, ySize=ySize, doApplyCalibration=True)
        tmplt_fake_img, tmplt_fake_sources = makeTestImage(
            psfSize=2.4, xSize=xSize, ySize=ySize, seed=9, nSrc=3, noiseLevel=0.25, fluxRange=1
        )
        # created fakes and added them to the images
        science.image += science_fake_img.image
        template.image += tmplt_fake_img.image

        # TODO: DM-40796 update to INJECTED names when source injection gets refactored
        # adding mask planes to both science and template images
        science.mask.addMaskPlane("FAKE")
        science_fake_bitmask = science.mask.getPlaneBitMask("FAKE")
        template.mask.addMaskPlane("FAKE")
        template_fake_bitmask = template.mask.getPlaneBitMask("FAKE")

        # makeTestImage sets the DETECTED plane on the sources; we can use
        # that to set the FAKE plane on the science and template images.
        detected = science_fake_img.mask.getPlaneBitMask("DETECTED")
        fake_pixels = (science_fake_img.mask.array & detected).nonzero()
        science.mask.array[fake_pixels] |= science_fake_bitmask
        detected = tmplt_fake_img.mask.getPlaneBitMask("DETECTED")
        fake_pixels = (tmplt_fake_img.mask.array & detected).nonzero()
        template.mask.array[fake_pixels] |= science_fake_bitmask

        science_fake_masked = (science.mask.array & science_fake_bitmask) > 0
        template_fake_masked = (template.mask.array & template_fake_bitmask) > 0

        subtractConfig = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        subtractConfig.sourceSelector.signalToNoise.fluxField = "truth_instFlux"
        subtractConfig.sourceSelector.signalToNoise.errField = "truth_instFluxErr"
        subtractTask = subtractImages.AlardLuptonSubtractTask(config=subtractConfig)
        subtraction = subtractTask.run(template, science, sources)

        # check subtraction mask plane is set where we set the previous masks
        diff_mask = subtraction.difference.mask

        # science mask should be now in INJECTED
        inj_masked = (diff_mask.array & diff_mask.getPlaneBitMask("INJECTED")) > 0

        # template mask should be now in INJECTED_TEMPLATE
        injTmplt_masked = (diff_mask.array & diff_mask.getPlaneBitMask("INJECTED_TEMPLATE")) > 0

        self.assertFloatsEqual(inj_masked.astype(int), science_fake_masked.astype(int))
        # The template is convolved, so the INJECTED_TEMPLATE mask plane may
        # include more pixels than the FAKE mask plane
        injTmplt_masked &= template_fake_masked
        self.assertFloatsEqual(injTmplt_masked.astype(int), template_fake_masked.astype(int))

        # Now check that detection of fakes have the correct flag for injections
        detectionTask = self._setup_detection()

        output = detectionTask.run(subtraction.matchedScience,
                                   subtraction.matchedTemplate,
                                   subtraction.difference,
                                   subtraction.kernelSources)
        diaSources = output.diaSources[~output.diaSources["sky_source"]].copy(deep=True)

        sci_refIds = []
        tmpl_refIds = []
        for diaSrc in diaSources:
            if diaSrc['base_PsfFlux_instFlux'] > 0:
                self._check_diaSource(science_fake_sources, diaSrc, scale=1, refIds=sci_refIds)
                self.assertTrue(diaSrc['base_PixelFlags_flag_injected'])
                self.assertTrue(diaSrc['base_PixelFlags_flag_injectedCenter'])
                self.assertFalse(diaSrc['base_PixelFlags_flag_injected_template'])
                self.assertFalse(diaSrc['base_PixelFlags_flag_injected_templateCenter'])
            else:
                self._check_diaSource(tmplt_fake_sources, diaSrc, scale=-1, refIds=tmpl_refIds)
                self.assertTrue(diaSrc['base_PixelFlags_flag_injected_template'])
                self.assertTrue(diaSrc['base_PixelFlags_flag_injected_templateCenter'])
                self.assertFalse(diaSrc['base_PixelFlags_flag_injected'])
                self.assertFalse(diaSrc['base_PixelFlags_flag_injectedCenter'])

    def test_mask_streaks(self):
        """Run detection on a difference image containing a streak.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 500
        xSize = 400
        ySize = 400
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "xSize": xSize, "ySize": ySize}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)

        # Configure the detection Task
        detectionTask = self._setup_detection(doMerge=False, doMaskStreaks=True,
                                              raiseOnNoDiaSources=False, raiseOnBadSubtractionRatio=False)

        # Test that no streaks are detected
        difference = science.clone()
        difference.maskedImage -= matchedTemplate.maskedImage
        output = detectionTask.run(science, matchedTemplate, difference, sources)
        outMask = output.subtractedMeasuredExposure.mask.array
        streakMask = output.subtractedMeasuredExposure.mask.getPlaneBitMask("STREAK")
        streakMaskSet = (outMask & streakMask) > 0
        self.assertTrue(np.all(streakMaskSet == 0))

        # Add streak-like shape and check that streak is detected
        difference.image.array[20:23, 40:200] += 50
        output = detectionTask.run(science, matchedTemplate, difference, sources)
        outMask = output.subtractedMeasuredExposure.mask.array
        streakMask = output.subtractedMeasuredExposure.mask.getPlaneBitMask("STREAK")
        streakMaskSet = (outMask & streakMask) > 0
        self.assertTrue(np.all(streakMaskSet[20:23, 40:200]))

        # Add a more trapezoid shaped streak across an image that is
        # fainter and check that it is also detected
        bbox = science.getBBox()
        difference = science.clone()
        difference.maskedImage -= matchedTemplate.maskedImage
        width = 100
        amplitude = 5
        x0 = -100 + bbox.getBeginX()
        y0 = -100 + bbox.getBeginY()
        x1 = xSize + x0 + 100
        y1 = ySize/2
        corner_coords = [lsst.geom.Point2D(x0, y0),
                         lsst.geom.Point2D(x0, y0 + width),
                         lsst.geom.Point2D(x1, y1),
                         lsst.geom.Point2D(x1 + width, y1)]
        streak_trapezoid = afwGeom.Polygon(corner_coords)
        streak_image = streak_trapezoid.createImage(bbox)
        streak_image *= amplitude
        difference.image.array += streak_image.array
        output = detectionTask.run(science, matchedTemplate, difference, sources)
        outMask = output.subtractedMeasuredExposure.mask.array
        streakMask = output.subtractedMeasuredExposure.mask.getPlaneBitMask("STREAK")
        streakMaskSet = (outMask & streakMask) > 0
        # Check that pixel values in streak_image equal the streak amplitude
        streak_check = np.where(streak_image.array == amplitude)
        self.assertTrue(np.all(streakMaskSet[streak_check]))
        # Check that the entire image was not masked STREAK
        self.assertFalse(np.all(streakMaskSet))

    def _setup_sattle_tests(self):
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 500
        shared_kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel,
                         "x0": 12345, "y0": 67890}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6,
                                         **shared_kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo(id=2))
        detector = DetectorWrapper(numAmps=1).detector
        science.setDetector(detector)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel / 4,
                                           noiseSeed=7, **shared_kwargs)
        difference = science.clone()

        detectionTask = self._setup_detection(doDeblend=True,
                                              badSubtractionRatioThreshold=1.,
                                              doSkySources=False, run_sattle=True)

        return science, matchedTemplate, difference, sources, detectionTask

    @mock.patch.dict(os.environ, {"SATTLE_URI_BASE": "fake_host:1234"})
    def test_sattle_not_available(self):
        science, matchedTemplate, difference, sources, detectionTask = self._setup_sattle_tests()

        response = MockResponse({"allow_list": []}, 500, "sattle internal error")
        with mock.patch('requests.put', return_value=response):
            with self.assertRaises(requests.exceptions.HTTPError):
                detectionTask.run(science, matchedTemplate, difference, sources,
                                  idFactory=IdFactory.makeSimple())

    @mock.patch.dict(os.environ, {"SATTLE_URI_BASE": "fake_host:1234"})
    def test_visit_id_not_in_sattle(self):
        science, matchedTemplate, difference, sources, detectionTask = self._setup_sattle_tests()

        response = MockResponse({"allow_list": []}, 404, "missing visit cache")
        # visit id not in sattle raises
        with self.assertRaises(requests.exceptions.HTTPError):
            with mock.patch('lsst.ip.diffim.detectAndMeasure.requests.put',
                            return_value=response):
                with mock.patch('lsst.ip.diffim.utils.populate_sattle_visit_cache'):
                    detectionTask.run(science, matchedTemplate, difference, sources,
                                      idFactory=IdFactory.makeSimple())

    @mock.patch.dict(os.environ, {"SATTLE_URI_BASE": "fake_host:1234"})
    def test_filter_satellites_some_allowed(self):
        science, matchedTemplate, difference, sources, detectionTask = self._setup_sattle_tests()

        allowed_ids = [1, 5]
        response = MockResponse({"allow_list": allowed_ids}, 200, "some allowed")
        with mock.patch('requests.put', return_value=response):
            output = detectionTask.run(science, matchedTemplate, difference, sources,
                                       idFactory=IdFactory.makeSimple())

        self.assertEqual(len(output.diaSources), 2)

        # Output should be sources 1 and 5 allowed out of 20
        self.assertEqual(set(output.diaSources['id']), set(allowed_ids))

    @mock.patch.dict(os.environ, {"SATTLE_URI_BASE": "fake_host:1234"})
    def test_filter_satellites_all_allowed(self):
        science, matchedTemplate, difference, sources, detectionTask = self._setup_sattle_tests()

        allowed_ids = list(range(1, 21))
        response = MockResponse({"allow_list": allowed_ids}, 200, "all allowed")
        # Run detection and check the results
        with mock.patch('requests.put', return_value=response):
            output = detectionTask.run(science, matchedTemplate, difference, sources,
                                       idFactory=IdFactory.makeSimple())

        ## Output should be all sources that went in. 20 go in, 20 should come out
        self.assertEqual(len(output.diaSources), 20)

        self.assertEqual(set(output.diaSources['id']), set(allowed_ids))

    @mock.patch.dict(os.environ, {"SATTLE_URI_BASE": "fake_host:1234"})
    def test_filter_satellites_none_allowed(self):
        science, matchedTemplate, difference, sources, detectionTask = self._setup_sattle_tests()

        response = MockResponse({"allow_list": []}, 200, "none allowed")
        # Run detection and confirm it raises for no diasources
        with self.assertRaises(detectAndMeasure.NoDiaSourcesError):
            with mock.patch('requests.put', return_value=response):
                detectionTask.run(science, matchedTemplate, difference, sources,
                                  idFactory=IdFactory.makeSimple())

    @mock.patch.dict(os.environ, {"SATTLE_URI_BASE": ""})
    def test_fail_on_sattle_misconfiguration(self):
        with self.assertRaises(pexConfig.FieldValidationError):
            self._setup_detection(run_sattle=True)

    def test_trailed_glints(self):
        """Test that the glint_trail column works, and that
        the trailed_glints output contains the expected information.
        """
        noiseLevel = 1.
        staticSeed = 1
        diffim, diaSources = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel, noiseSeed=6)
        self._check_values(diaSources['glint_trail'])

        # Run detection and return the output Struct so we can check it
        def _detection_wrapper(diffim, diaSources):
            detectionTask = self._setup_detection()
            scienceBase, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6)
            matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7)
            science = scienceBase.clone()
            science.maskedImage -= diffim.maskedImage
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            output = detectionTask.run(science, matchedTemplate, difference, sources)
            return output

        output = _detection_wrapper(diffim, diaSources)
        self.assertTrue('slopes' in output.glintTrailInfo)
        self.assertTrue('intercepts' in output.glintTrailInfo)
        self.assertTrue('stderrs' in output.glintTrailInfo)
        self.assertTrue('lengths' in output.glintTrailInfo)
        self.assertTrue('angles' in output.glintTrailInfo)


class DetectAndMeasureScoreTest(DetectAndMeasureTestBase, lsst.utils.tests.TestCase):
    detectionTask = detectAndMeasure.DetectAndMeasureScoreTask

    def test_detection_xy0(self):
        """Basic functionality test with non-zero x0 and y0.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "x0": 12345, "y0": 67890}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()
        subtractTask = subtractImages.AlardLuptonPreconvolveSubtractTask()
        scienceKernel = science.psf.getKernel()
        score = subtractTask._convolveExposure(difference, scienceKernel, subtractTask.convolutionControl)

        # Configure the detection Task
        # Set the subtraction residual metric threshold high, since the template
        # is not actually subtracted from the science image.
        detectionTask = self._setup_detection(doDeblend=True, badSubtractionRatioThreshold=1.,
                                              doSkySources=False)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference, score, sources,
                                   idFactory=self.idGenerator.make_table_id_factory())

        # Catalog ids should be very large from this id generator.
        self.assertTrue(all(output.diaSources['id'] > 1000000000))
        subtractedMeasuredExposure = output.subtractedMeasuredExposure

        self.assertImagesEqual(subtractedMeasuredExposure.image, difference.image)

        # Not all of the sources will be detected: preconvolution results in
        # a larger edge mask, so we miss an edge source.
        self.assertEqual(len(output.diaSources), len(sources)-1)
        # no sources should be flagged as negative
        self.assertEqual(len(~output.diaSources["is_negative"]), len(output.diaSources))
        # TODO DM-41496: restore this block once we handle detections on edge
        # pixels better; at least one of these sources currently has a bad
        # centroid because most of the source is rejected as EDGE.
        # refIds = []
        # for source in sources:
        #     self._check_diaSource(output.diaSources, source, refIds=refIds)

    def test_measurements_finite(self):
        """Measured fluxes and centroids should always be finite.
        """
        columnNames = ["coord_ra", "coord_dec", "ip_diffim_forced_PsfFlux_instFlux"]

        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        xSize = 256
        ySize = 256
        kwargs = {"psfSize": 2.4, "x0": 0, "y0": 0,
                  "xSize": xSize, "ySize": ySize}
        science, sources = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel, noiseSeed=6,
                                         nSrc=1, **kwargs)
        matchedTemplate, _ = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel/4, noiseSeed=7,
                                           nSrc=1, **kwargs)
        rng = np.random.RandomState(3)
        xLoc = np.arange(-5, xSize+5, 10)
        rng.shuffle(xLoc)
        yLoc = np.arange(-5, ySize+5, 10)
        rng.shuffle(yLoc)
        transients, transientSources = makeTestImage(seed=transientSeed,
                                                     nSrc=len(xLoc), fluxLevel=1000.,
                                                     noiseLevel=noiseLevel, noiseSeed=8,
                                                     xLoc=xLoc, yLoc=yLoc,
                                                     **kwargs)
        difference = science.clone()
        difference.maskedImage -= matchedTemplate.maskedImage
        difference.maskedImage += transients.maskedImage
        subtractTask = subtractImages.AlardLuptonPreconvolveSubtractTask()
        scienceKernel = science.psf.getKernel()
        score = subtractTask._convolveExposure(difference, scienceKernel, subtractTask.convolutionControl)

        # Configure the detection Task
        detectionTask = self._setup_detection(doForcedMeasurement=True)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference, score, sources)

        for column in columnNames:
            self._check_values(output.diaSources[column])
        self._check_values(output.diaSources.getX(), minValue=0, maxValue=xSize)
        self._check_values(output.diaSources.getY(), minValue=0, maxValue=ySize)
        self._check_values(output.diaSources.getPsfInstFlux())

    def test_detect_transients(self):
        """Run detection on a difference image containing transients.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        scienceBase, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        scienceKernel = scienceBase.psf.getKernel()
        subtractTask = subtractImages.AlardLuptonPreconvolveSubtractTask()

        # Configure the detection Task
        detectionTask = self._setup_detection(doMerge=False)
        kwargs["seed"] = transientSeed
        kwargs["nSrc"] = 10
        kwargs["fluxLevel"] = 1000

        # Run detection and check the results
        def _detection_wrapper(positive=True):
            """Simulate positive or negative transients and run detection.

            Parameters
            ----------
            positive : `bool`, optional
                If set, use positive transient sources.
            """

            transients, transientSources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=8, **kwargs)
            science = scienceBase.clone()
            if positive:
                science.maskedImage += transients.maskedImage
            else:
                science.maskedImage -= transients.maskedImage
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            score = subtractTask._convolveExposure(difference, scienceKernel, subtractTask.convolutionControl)
            # NOTE: NoiseReplacer (run by forcedMeasurement) can modify the
            # science image if we've e.g. removed parents post-deblending.
            # Pass a clone of the science image, so that it doesn't disrupt
            # later tests.
            output = detectionTask.run(science, matchedTemplate, difference, score, sources)
            refIds = []
            scale = 1. if positive else -1.
            # sources near the edge may have untrustworthy centroids
            goodSrcFlags = ~output.diaSources['base_PixelFlags_flag_edge']
            for diaSource, goodSrcFlag in zip(output.diaSources, goodSrcFlags):
                if goodSrcFlag:
                    self._check_diaSource(transientSources, diaSource, refIds=refIds, scale=scale)
        _detection_wrapper(positive=True)
        _detection_wrapper(positive=False)

    def test_detect_dipoles(self):
        """Run detection on a difference image containing dipoles.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 1000
        fluxRange = 1.5
        nSources = 10
        offset = 1
        xSize = 300
        ySize = 300
        kernelSize = 32
        # Avoid placing sources near the edge for this test, so that we can
        # easily check that the correct number of sources are detected.
        templateBorderSize = kernelSize//2
        dipoleFlag = "ip_diffim_DipoleFit_classification"
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "fluxRange": fluxRange,
                  "nSrc": nSources, "templateBorderSize": templateBorderSize, "kernelSize": kernelSize,
                  "xSize": xSize, "ySize": ySize}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()
        # Shift the template by a pixel in order to make dipoles in the difference image.
        matchedTemplate.image.array[...] = np.roll(matchedTemplate.image.array[...], offset, axis=0)
        matchedTemplate.variance.array[...] = np.roll(matchedTemplate.variance.array[...], offset, axis=0)
        matchedTemplate.mask.array[...] = np.roll(matchedTemplate.mask.array[...], offset, axis=0)
        difference.maskedImage -= matchedTemplate.maskedImage[science.getBBox()]
        subtractTask = subtractImages.AlardLuptonPreconvolveSubtractTask()
        scienceKernel = science.psf.getKernel()
        score = subtractTask._convolveExposure(difference, scienceKernel, subtractTask.convolutionControl)

        detectionTask = self._setup_detection(badSubtractionRatioThreshold=0.3)
        output = detectionTask.run(science, matchedTemplate, difference, score, sources)
        diaSources = output.diaSources[~output.diaSources["sky_source"]].copy(deep=True)
        self.assertEqual(len(diaSources), len(sources))
        refIds = []
        for diaSource in diaSources:
            if diaSource[dipoleFlag]:
                self._check_diaSource(sources, diaSource, refIds=refIds, scale=0,
                                      rtol=0.05, atol=None, usePsfFlux=False)
                self.assertFloatsAlmostEqual(diaSource["ip_diffim_DipoleFit_orientation"],
                                             -np.pi / 2, atol=2.)
                self.assertFloatsAlmostEqual(diaSource["ip_diffim_DipoleFit_separation"], offset, rtol=0.1)
            else:
                raise ValueError("DiaSource with ID %s is not a dipole!", diaSource.getId())

    def test_sky_sources(self):
        """Add sky sources and check that they are sufficiently far from other
        sources and have negligible flux.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        transientFluxLevel = 1000.
        transientFluxRange = 1.5
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        transients, transientSources = makeTestImage(seed=transientSeed, psfSize=2.4,
                                                     nSrc=10, fluxLevel=transientFluxLevel,
                                                     fluxRange=transientFluxRange,
                                                     noiseLevel=noiseLevel, noiseSeed=8)
        difference = science.clone()
        difference.maskedImage -= matchedTemplate.maskedImage
        difference.maskedImage += transients.maskedImage
        subtractTask = subtractImages.AlardLuptonPreconvolveSubtractTask()
        scienceKernel = science.psf.getKernel()
        kernelWidth = np.max(scienceKernel.getDimensions())//2
        score = subtractTask._convolveExposure(difference, scienceKernel, subtractTask.convolutionControl)

        # Configure the detection Task
        detectionTask = self._setup_detection(doSkySources=True)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference, score, sources,
                                   idFactory=self.idGenerator.make_table_id_factory())
        nSkySourcesGenerated = detectionTask.metadata["n_skySources"]
        skySources = output.diaSources[output.diaSources["sky_source"]]
        self.assertEqual(len(skySources), nSkySourcesGenerated)
        for skySource in skySources:
            # Disable the sky_source flag to enable checks.
            skySource["sky_source"] = False
            # The sky sources should not be close to any other source
            with self.assertRaises(AssertionError):
                self._check_diaSource(transientSources, skySource, matchDistance=kernelWidth)
            with self.assertRaises(AssertionError):
                self._check_diaSource(sources, skySource, matchDistance=kernelWidth)
            # The sky sources should have low flux levels.
            self._check_diaSource(transientSources, skySource, matchDistance=1000, scale=0.,
                                  atol=np.sqrt(transientFluxRange*transientFluxLevel))

        # Catalog ids should be very large from this id generator.
        self.assertTrue(all(output.diaSources['id'] > 1000000000))

    def test_exclude_mask_detections(self):
        """Sources with certain bad mask planes set should not be detected.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        fluxLevel = 500
        radius = 2
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        science.getInfo().setVisitInfo(makeVisitInfo())
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)

        subtractTask = subtractImages.AlardLuptonPreconvolveSubtractTask()
        scienceKernel = science.psf.getKernel()
        # Configure the detection Task
        detectionTask = self._setup_detection()
        excludeMaskPlanes = detectionTask.config.detection.excludeMaskPlanes
        nBad = len(excludeMaskPlanes)
        self.assertGreater(nBad, 0)
        kwargs["seed"] = transientSeed
        kwargs["nSrc"] = nBad
        kwargs["fluxLevel"] = 1000

        # Run detection and check the results
        def _detection_wrapper(setFlags=True):
            transients, transientSources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=8, **kwargs)
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            difference.maskedImage += transients.maskedImage
            if setFlags:
                for src, badMask in zip(transientSources, excludeMaskPlanes):
                    srcX = int(src.getX())
                    srcY = int(src.getY())
                    srcBbox = lsst.geom.Box2I(lsst.geom.Point2I(srcX - radius, srcY - radius),
                                              lsst.geom.Extent2I(2*radius + 1, 2*radius + 1))
                    difference[srcBbox].mask.array |= lsst.afw.image.Mask.getPlaneBitMask(badMask)
            score = subtractTask._convolveExposure(difference, scienceKernel, subtractTask.convolutionControl)
            if setFlags:
                with self.assertRaises(AlgorithmError):
                    output = detectionTask.run(science, matchedTemplate, difference, score, sources)
                return
            else:
                output = detectionTask.run(science, matchedTemplate, difference, score, sources)
                refIds = []
                goodSrcFlags = checkMask(difference.mask, transientSources, excludeMaskPlanes)
                self.assertEqual(np.sum(~goodSrcFlags), 0)
                for diaSource, goodSrcFlag in zip(output.diaSources, goodSrcFlags):
                    if ~goodSrcFlag:
                        with self.assertRaises(AssertionError):
                            self._check_diaSource(transientSources, diaSource, refIds=refIds)
                    else:
                        self._check_diaSource(transientSources, diaSource, refIds=refIds)
        _detection_wrapper(setFlags=False)
        _detection_wrapper(setFlags=True)


class TestNegativePeaks(lsst.utils.tests.TestCase):
    """Tests of deblending and merging negative peaks, to test fixes for the
    various problems found on DM-48596.
    """

    def testDeblendNegatives(self):
        """Test that negative peaks get deblended and not destroyed: DM-48704.
        This is only a test of deblending, not of merging.
        """
        # Make a science image with one blend of two positive sources, to
        # subtract from an empty template, resulting in a negative diffim.
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(100, 100))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        delta = 10
        with dataset.addBlend() as family:
            family.addChild(instFlux=2E5, centroid=lsst.geom.Point2D(50, 72))
            family.addChild(instFlux=2.5E5, centroid=lsst.geom.Point2D(50+delta, 74))
        science, catalog = dataset.realize(noise=1.0,
                                           schema=lsst.meas.base.tests.TestDataset.makeMinimalSchema())
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        template, _ = dataset.realize(noise=1.0,
                                      schema=lsst.meas.base.tests.TestDataset.makeMinimalSchema())
        difference = template.clone()
        difference.image -= science.image

        config = detectAndMeasure.DetectAndMeasureTask.ConfigClass()
        config.doDeblend = True
        task = detectAndMeasure.DetectAndMeasureTask(config=config)
        # prelude steps taken from `detectAndMeasure.run`
        task._prepareInputs(difference)
        table = lsst.afw.table.SourceTable.make(task.schema)
        results = task.detection.run(table=table, exposure=difference, doSmooth=True)
        # Just run the deblend step so we can check the footprints independently.
        sources, positives, negatives = task._deblend(difference, results.positive, results.negative)

        # DM-48704 fixed a problem where the peaks were in the footprints, but
        # the spans were empty.
        self.assertEqual(len(negatives), 2)
        self.assertEqual(len(positives), 0)
        self.assertGreater(negatives[0].getFootprint().getSpans().getArea(), 0)
        self.assertGreater(negatives[1].getFootprint().getSpans().getArea(), 0)
        # Deblended children are HeavyFootprints; we have to make sure the
        # pixel values in those are correct (though DetectAndMeasureTask
        # doesn't use the fact that they're Heavy).
        # The sources are positive in the science image, and negative in the
        # diffim, so the minimum value in the deblended negative footprint is
        # the maximum value in the science catalog footprint (ignoring the
        # noise in the science and template, hence rtol).
        # (catalog[0] is the parent; we want the children)
        self.assertFloatsAlmostEqual(negatives[0].getFootprint().getImageArray().min(),
                                     -catalog[1].getFootprint().getImageArray().max(), rtol=1e-4)
        self.assertFloatsAlmostEqual(negatives[1].getFootprint().getImageArray().min(),
                                     -catalog[2].getFootprint().getImageArray().max(), rtol=1e-4)

        self.assertEqual(sources["is_negative"].sum(), 2)

    def testMergeFootprints(self):
        """Test that merging footprints does not cause negative ones to
        disappear (e.g. get merged into non-connected footprints).

        As implemented, the diffim will have three positive sources (one a
        blended pair), and 7 negative sources (two a blended pair).
        """
        # Make a template image with multiple blends, designed to trigger the
        # negative-footprint-related bug in `FootprintSet.merge`.
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(400, 200))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        # Because these sources are on the template, they will be negative on
        # the diffim, unless the pixels are explicitly made negative via
        # `template.image.subset` below.
        # positive isolated source on diffim
        dataset.addSource(instFlux=.5E5, centroid=lsst.geom.Point2D(25, 26))
        # negative isolated source on diffim
        dataset.addSource(instFlux=.7E5, centroid=lsst.geom.Point2D(75, 24),
                          shape=lsst.afw.geom.Quadrupole(12, 7, 2))
        delta = 10
        # negative blended pair on diffim
        with dataset.addBlend() as family:
            family.addChild(instFlux=1E5, centroid=lsst.geom.Point2D(150, 72))
            family.addChild(instFlux=1.5E5, centroid=lsst.geom.Point2D(150+delta, 74))
        # positive blended pair on diffim
        with dataset.addBlend() as family:
            family.addChild(instFlux=2E5, centroid=lsst.geom.Point2D(250, 72))
            family.addChild(instFlux=2.5E5, centroid=lsst.geom.Point2D(250+delta, 74))
        # negative blended pair on diffim
        with dataset.addBlend() as family:
            family.addChild(instFlux=3E5, centroid=lsst.geom.Point2D(350, 72))
            family.addChild(instFlux=3.5E5, centroid=lsst.geom.Point2D(350+delta, 74))
        # negative isolated source on diffim
        dataset.addSource(instFlux=4E5, centroid=lsst.geom.Point2D(75, 124))
        # negative isolated source on diffim
        dataset.addSource(instFlux=5E5, centroid=lsst.geom.Point2D(175, 124))
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        schema.addField("sky_source", "Flag", "Sky source.")
        template, catalog = dataset.realize(noise=100.0, schema=schema)
        # These will be positive sources on the diffim (as noted above).
        template.image.subset(lsst.geom.Box2I(lsst.geom.Point2I(15, 15),
                                              lsst.geom.Point2I(35, 35))).array *= -1
        template.image.subset(lsst.geom.Box2I(lsst.geom.Point2I(230, 60),
                                              lsst.geom.Point2I(270, 85))).array *= -1
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        science, _ = dataset.realize(noise=100.0, schema=schema)
        difference = science.clone()
        difference.image -= template.image

        config = detectAndMeasure.DetectAndMeasureTask.ConfigClass()
        config.doDeblend = True
        config.raiseOnBadSubtractionRatio = False
        config.raiseOnNoDiaSources = False
        task = detectAndMeasure.DetectAndMeasureTask(config=config)
        result = task.run(science, template, difference, catalog)

        # The original bug manifested as the (175,124) source disappaering in
        # the merge step, and being included in the peaks of a duplicated
        # (75,124) source, even though their footprints are not connected.
        mask = np.isclose(result.diaSources["slot_Centroid_x"], 75) & \
            np.isclose(result.diaSources["slot_Centroid_y"], 124)
        self.assertEqual(mask.sum(), 1)
        peaks = result.diaSources[mask][0].getFootprint().peaks
        self.assertEqual(len(peaks), 1)
        self.assertEqual(peaks[0].getIx(), 75)
        self.assertEqual(peaks[0].getIy(), 124)

        mask = np.isclose(result.diaSources["slot_Centroid_x"], 175) & \
            np.isclose(result.diaSources["slot_Centroid_y"], 124)
        self.assertEqual(mask.sum(), 1)
        peaks = result.diaSources[mask][0].getFootprint().peaks
        self.assertEqual(len(peaks), 1)
        self.assertEqual(peaks[0].getIx(), 175)
        self.assertEqual(peaks[0].getIy(), 124)

        # This checks that all the returned footprints have exactly one peak;
        # it's a more general test than the above, but I'm keeping that for
        # the more explicit test of the original bugged sources.
        peak_x = np.array([peak['f_x'] for src in result.diaSources for peak in src.getFootprint().peaks])
        peak_y = np.array([peak['f_y'] for src in result.diaSources for peak in src.getFootprint().peaks])
        peaks = np.column_stack((peak_x, peak_y))
        unique_peaks, counts = np.unique(peaks, axis=0, return_counts=True)
        self.assertEqual(np.sum(counts > 1), 0)

        # There are three positive sources that should not have `is_negative`
        # set, independent of how deblending/merging of negative sources is
        # handled.
        self.assertEqual((~result.diaSources["is_negative"]).sum(), 3)


def makeVisitInfo(id=1):
    """Return a non-NaN visitInfo."""
    return afwImage.VisitInfo(id=id,
                              exposureTime=10.01,
                              darkTime=11.02,
                              date=dafBase.DateTime(65321.1, dafBase.DateTime.MJD, dafBase.DateTime.TAI),
                              ut1=12345.1,
                              era=45.1*geom.degrees,
                              boresightRaDec=geom.SpherePoint(23.1, 73.2, geom.degrees),
                              boresightAzAlt=geom.SpherePoint(134.5, 33.3, geom.degrees),
                              boresightAirmass=1.73,
                              boresightRotAngle=73.2*geom.degrees,
                              rotType=afwImage.RotType.SKY,
                              observatory=Observatory(
                                  11.1*geom.degrees, 22.2*geom.degrees, 0.333),
                              weather=Weather(1.1, 2.2, 34.5),
                              )


class MockResponse:
    """Provide a mock for requests.put calls"""
    def __init__(self, json_data, status_code, text):
        self.json_data = json_data
        self.status_code = status_code
        self.text = text

    def json(self):
        return self.json_data

    def raise_for_status(self):
        if self.status_code != 200:
            raise requests.exceptions.HTTPError


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
