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
import unittest

import lsst.geom
from lsst.ip.diffim import detectAndMeasure, subtractImages
from lsst.ip.diffim.utils import makeTestImage
import lsst.utils.tests


class DetectAndMeasureTestBase(lsst.utils.tests.TestCase):

    def _check_diaSource(self, refSources, diaSource, refIds=None,
                         matchDistance=1., scale=1., usePsfFlux=True,
                         rtol=0.02, atol=None):
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

    def _setup_detection(self, doSkySources=False, nSkySources=5, **kwargs):
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
        return self.detectionTask(config=config)


class DetectAndMeasureTest(DetectAndMeasureTestBase):
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
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()

        # Configure the detection Task
        detectionTask = self._setup_detection()

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference)
        subtractedMeasuredExposure = output.subtractedMeasuredExposure

        self.assertImagesEqual(subtractedMeasuredExposure.image, difference.image)

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

        # Configure the detection Task
        detectionTask = self._setup_detection(doForcedMeasurement=True)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference)

        for column in columnNames:
            self._check_values(output.diaSources[column])
        self._check_values(output.diaSources.getX(), minValue=0, maxValue=xSize)
        self._check_values(output.diaSources.getY(), minValue=0, maxValue=ySize)
        self._check_values(output.diaSources.getPsfInstFlux())

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
        matchedTemplate, _ = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel/4, noiseSeed=7,
                                           nSrc=1, **kwargs)
        difference = science.clone()
        bbox = difference.getBBox()
        difference.maskedImage -= matchedTemplate.maskedImage

        # Configure the detection Task, and do not remove unphysical sources
        detectionTask = self._setup_detection(doForcedMeasurement=False, doSkySources=True, nSkySources=20,
                                              badSourceFlags=["base_PixelFlags_flag_offimage", ])

        # Run detection and check the results
        diaSources = detectionTask.run(science, matchedTemplate, difference).diaSources
        badDiaSrc0 = ~bbox.contains(diaSources.getX(), diaSources.getY())
        nBad0 = np.count_nonzero(badDiaSrc0)
        # Verify that all sources are physical
        self.assertEqual(nBad0, 0)
        # Set a few centroids outside the image bounding box
        nSetBad = 5
        for src in diaSources[0: nSetBad]:
            src["slot_Centroid_x"] += xSize
            src["slot_Centroid_y"] += ySize
            src["base_PixelFlags_flag_offimage"] = True
        # Verify that these sources are outside the image
        badDiaSrc1 = ~bbox.contains(diaSources.getX(), diaSources.getY())
        nBad1 = np.count_nonzero(badDiaSrc1)
        self.assertEqual(nBad1, nSetBad)
        diaSources2 = detectionTask.removeBadSources(diaSources)
        badDiaSrc2 = ~bbox.contains(diaSources2.getX(), diaSources2.getY())
        nBad2 = np.count_nonzero(badDiaSrc2)

        # Verify that no sources outside the image bounding box remain
        self.assertEqual(nBad2, 0)
        self.assertEqual(len(diaSources2), len(diaSources) - nSetBad)

    def test_remove_nan_centroid(self):
        """Check that sources with non-finite centroids are removed from the catalog.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        xSize = 256
        ySize = 256
        kwargs = {"psfSize": 2.4, "xSize": xSize, "ySize": ySize}
        science, sources = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel, noiseSeed=6,
                                         nSrc=1, **kwargs)
        matchedTemplate, _ = makeTestImage(seed=staticSeed, noiseLevel=noiseLevel/4, noiseSeed=7,
                                           nSrc=1, **kwargs)
        difference = science.clone()
        bbox = difference.getBBox()
        difference.maskedImage -= matchedTemplate.maskedImage

        # Configure the detection Task, and do not remove unphysical sources
        detectionTask = self._setup_detection(doForcedMeasurement=False, doSkySources=True, nSkySources=20,
                                              badSourceFlags=["base_PixelFlags_flag_offimage", ])

        # Run detection and check the results
        diaSources = detectionTask.run(science, matchedTemplate, difference).diaSources
        badDiaSrc0 = ~bbox.contains(diaSources.getX(), diaSources.getY())
        nBad0 = np.count_nonzero(badDiaSrc0)
        # Verify that all sources are physical
        self.assertEqual(nBad0, 0)
        # Set a few centroids outside the image bounding box
        nSetBad = 5
        for i, src in enumerate(diaSources[0: nSetBad]):
            if i % 3 == 0:
                src["slot_Centroid_x"] = np.nan
            elif i % 3 == 1:
                src["slot_Centroid_y"] = np.nan
            elif i % 3 == 2:
                src["slot_Centroid_x"] = np.nan
                src["slot_Centroid_y"] = np.nan
        # Verify that these sources are outside the image
        badDiaSrc1 = ~bbox.contains(diaSources.getX(), diaSources.getY())
        nBad1 = np.count_nonzero(badDiaSrc1)
        self.assertEqual(nBad1, nSetBad)
        diaSources2 = detectionTask.removeBadSources(diaSources)
        badDiaSrc2 = ~bbox.contains(diaSources2.getX(), diaSources2.getY())
        nBad2 = np.count_nonzero(badDiaSrc2)

        # Verify that no sources outside the image bounding box remain
        self.assertEqual(nBad2, 0)
        self.assertEqual(len(diaSources2), len(diaSources) - nSetBad)

    def test_detect_transients(self):
        """Run detection on a difference image containing transients.
        """
        # Set up the simulated images
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)

        # Configure the detection Task
        detectionTask = self._setup_detection(doMerge=False)
        kwargs["seed"] = transientSeed
        kwargs["nSrc"] = 10
        kwargs["fluxLevel"] = 1000

        # Run detection and check the results
        def _detection_wrapper(positive=True):
            transients, transientSources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=8, **kwargs)
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            if positive:
                difference.maskedImage += transients.maskedImage
            else:
                difference.maskedImage -= transients.maskedImage
            output = detectionTask.run(science, matchedTemplate, difference)
            refIds = []
            scale = 1. if positive else -1.
            for diaSource in output.diaSources:
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
        dipoleFlag = "ip_diffim_DipoleFit_flag_classification"
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "fluxRange": fluxRange,
                  "nSrc": nSources, "templateBorderSize": templateBorderSize, "kernelSize": kernelSize,
                  "xSize": xSize, "ySize": ySize}
        dipoleFlag = "ip_diffim_DipoleFit_flag_classification"
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()
        matchedTemplate.image.array[...] = np.roll(matchedTemplate.image.array[...], offset, axis=0)
        matchedTemplate.variance.array[...] = np.roll(matchedTemplate.variance.array[...], offset, axis=0)
        matchedTemplate.mask.array[...] = np.roll(matchedTemplate.mask.array[...], offset, axis=0)
        difference.maskedImage -= matchedTemplate.maskedImage[science.getBBox()]

        # Configure the detection Task
        detectionTask = self._setup_detection(doMerge=False)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference)
        self.assertIn(dipoleFlag, output.diaSources.schema.getNames())
        nSourcesDet = len(sources)
        self.assertEqual(len(output.diaSources), 2*nSourcesDet)
        refIds = []
        # The diaSource check should fail if we don't merge positive and negative footprints
        for diaSource in output.diaSources:
            with self.assertRaises(AssertionError):
                self._check_diaSource(sources, diaSource, refIds=refIds, scale=0,
                                      atol=np.sqrt(fluxRange*fluxLevel))

        detectionTask2 = self._setup_detection(doMerge=True)
        output2 = detectionTask2.run(science, matchedTemplate, difference)
        self.assertEqual(len(output2.diaSources), nSourcesDet)
        refIds = []
        for diaSource in output2.diaSources:
            if diaSource[dipoleFlag]:
                self._check_diaSource(sources, diaSource, refIds=refIds, scale=0,
                                      rtol=0.05, atol=None, usePsfFlux=False)
                self.assertFloatsAlmostEqual(diaSource["ip_diffim_DipoleFit_orientation"], -90., atol=2.)
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
        kernelWidth = np.max(science.psf.getKernel().getDimensions())//2

        # Configure the detection Task
        detectionTask = self._setup_detection(doSkySources=True)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference)
        skySources = output.diaSources[output.diaSources["sky_source"]]
        self.assertEqual(len(skySources), detectionTask.config.skySources.nSources)
        for skySource in skySources:
            # The sky sources should not be close to any other source
            with self.assertRaises(AssertionError):
                self._check_diaSource(transientSources, skySource, matchDistance=kernelWidth)
            with self.assertRaises(AssertionError):
                self._check_diaSource(sources, skySource, matchDistance=kernelWidth)
            # The sky sources should have low flux levels.
            self._check_diaSource(transientSources, skySource, matchDistance=1000, scale=0.,
                                  atol=np.sqrt(transientFluxRange*transientFluxLevel))

    def test_edge_detections(self):
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
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)

        _checkMask = subtractImages.AlardLuptonSubtractTask._checkMask
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
            output = detectionTask.run(science, matchedTemplate, difference)
            refIds = []
            goodSrcFlags = _checkMask(difference.mask, transientSources, excludeMaskPlanes)
            if setFlags:
                self.assertEqual(np.sum(~goodSrcFlags), nBad)
            else:
                self.assertEqual(np.sum(~goodSrcFlags), 0)
            for diaSource, goodSrcFlag in zip(output.diaSources, goodSrcFlags):
                if ~goodSrcFlag:
                    with self.assertRaises(AssertionError):
                        self._check_diaSource(transientSources, diaSource, refIds=refIds)
                else:
                    self._check_diaSource(transientSources, diaSource, refIds=refIds)
        _detection_wrapper(setFlags=False)
        _detection_wrapper(setFlags=True)


class DetectAndMeasureScoreTest(DetectAndMeasureTestBase):
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
        detectionTask = self._setup_detection()

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference, score)
        subtractedMeasuredExposure = output.subtractedMeasuredExposure

        self.assertImagesEqual(subtractedMeasuredExposure.image, difference.image)

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
        output = detectionTask.run(science, matchedTemplate, difference, score)

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
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        scienceKernel = science.psf.getKernel()
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
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            if positive:
                difference.maskedImage += transients.maskedImage
            else:
                difference.maskedImage -= transients.maskedImage
            score = subtractTask._convolveExposure(difference, scienceKernel, subtractTask.convolutionControl)
            output = detectionTask.run(science, matchedTemplate, difference, score)
            refIds = []
            scale = 1. if positive else -1.
            goodSrcFlags = subtractTask._checkMask(score.mask, transientSources,
                                                   subtractTask.config.badMaskPlanes)
            for diaSource, goodSrcFlag in zip(output.diaSources, goodSrcFlags):
                if ~goodSrcFlag:
                    with self.assertRaises(AssertionError):
                        self._check_diaSource(transientSources, diaSource, refIds=refIds, scale=scale)
                else:
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
        dipoleFlag = "ip_diffim_DipoleFit_flag_classification"
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

        # Configure the detection Task
        detectionTask = self._setup_detection(doMerge=False)

        # Run detection and check the results
        output = detectionTask.run(science, matchedTemplate, difference, score)
        self.assertIn(dipoleFlag, output.diaSources.schema.getNames())
        nSourcesDet = len(sources)
        # Since we did not merge the dipoles, each source should result in
        # both a positive and a negative diaSource
        self.assertEqual(len(output.diaSources), 2*nSourcesDet)
        refIds = []
        # The diaSource check should fail if we don't merge positive and negative footprints
        for diaSource in output.diaSources:
            with self.assertRaises(AssertionError):
                self._check_diaSource(sources, diaSource, refIds=refIds, scale=0,
                                      atol=np.sqrt(fluxRange*fluxLevel))

        detectionTask2 = self._setup_detection(doMerge=True)
        output2 = detectionTask2.run(science, matchedTemplate, difference, score)
        self.assertEqual(len(output2.diaSources), nSourcesDet)
        refIds = []
        for diaSource in output2.diaSources:
            if diaSource[dipoleFlag]:
                self._check_diaSource(sources, diaSource, refIds=refIds, scale=0,
                                      rtol=0.05, atol=None, usePsfFlux=False)
                self.assertFloatsAlmostEqual(diaSource["ip_diffim_DipoleFit_orientation"], -90., atol=2.)
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
        output = detectionTask.run(science, matchedTemplate, difference, score)
        skySources = output.diaSources[output.diaSources["sky_source"]]
        self.assertEqual(len(skySources), detectionTask.config.skySources.nSources)
        for skySource in skySources:
            # The sky sources should not be close to any other source
            with self.assertRaises(AssertionError):
                self._check_diaSource(transientSources, skySource, matchDistance=kernelWidth)
            with self.assertRaises(AssertionError):
                self._check_diaSource(sources, skySource, matchDistance=kernelWidth)
            # The sky sources should have low flux levels.
            self._check_diaSource(transientSources, skySource, matchDistance=1000, scale=0.,
                                  atol=np.sqrt(transientFluxRange*transientFluxLevel))

    def test_edge_detections(self):
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
            output = detectionTask.run(science, matchedTemplate, difference, score)
            refIds = []
            goodSrcFlags = subtractTask._checkMask(difference.mask, transientSources, excludeMaskPlanes)
            if setFlags:
                self.assertEqual(np.sum(~goodSrcFlags), nBad)
            else:
                self.assertEqual(np.sum(~goodSrcFlags), 0)
            for diaSource, goodSrcFlag in zip(output.diaSources, goodSrcFlags):
                if ~goodSrcFlag:
                    with self.assertRaises(AssertionError):
                        self._check_diaSource(transientSources, diaSource, refIds=refIds)
                else:
                    self._check_diaSource(transientSources, diaSource, refIds=refIds)
        _detection_wrapper(setFlags=False)
        _detection_wrapper(setFlags=True)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
