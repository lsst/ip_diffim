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
from lsst.ip.diffim import detectAndMeasure
from lsst.ip.diffim.utils import makeTestImage
import lsst.utils.tests


class DetectAndMeasureTest(lsst.utils.tests.TestCase):

    def test_detection_runs(self):
        """Basic smoke test.
        """
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "x0": 0, "y0": 0}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()
        config = detectAndMeasure.DetectAndMeasureTask.ConfigClass()
        config.doApCorr = False
        task = detectAndMeasure.DetectAndMeasureTask(config=config)
        output = task.run(science, matchedTemplate, difference)
        subtractedMeasuredExposure = output.subtractedMeasuredExposure
        self.assertImagesEqual(subtractedMeasuredExposure.image, difference.image)

    def test_detection_xy0(self):
        """Basic smoke test with non-zero x0 and y0.
        """
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "x0": 12345, "y0": 67890}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()
        config = detectAndMeasure.DetectAndMeasureTask.ConfigClass()
        config.doApCorr = False
        config.doMerge = False
        config.doSkySources = False
        config.doForcedMeasurement = False
        task = detectAndMeasure.DetectAndMeasureTask(config=config)
        output = task.run(science, matchedTemplate, difference)
        subtractedMeasuredExposure = output.subtractedMeasuredExposure

        self.assertImagesEqual(subtractedMeasuredExposure.image, difference.image)

    def test_detect_transients(self):
        """Run detection on a difference image containing transients.
        """
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        config = detectAndMeasure.DetectAndMeasureTask.ConfigClass()
        config.doApCorr = False
        config.doMerge = False
        config.doSkySources = False
        config.doForcedMeasurement = False
        detectionTask = detectAndMeasure.DetectAndMeasureTask(config=config)

        def _detection_wrapper(polarity=1):
            transients, transientSources = makeTestImage(seed=transientSeed, psfSize=2.4,
                                                         nSrc=10, fluxLevel=1000.,
                                                         noiseLevel=noiseLevel, noiseSeed=8)
            difference = science.clone()
            difference.maskedImage -= matchedTemplate.maskedImage
            if polarity > 0:
                difference.maskedImage += transients.maskedImage
            else:
                difference.maskedImage -= transients.maskedImage
            output = detectionTask.run(science, matchedTemplate, difference)
            refIds = []
            for diaSource in output.diaSources:
                self._check_diaSource(transientSources, diaSource, refIds=refIds, scale=polarity)
        _detection_wrapper(polarity=1)
        _detection_wrapper(polarity=-1)

    def test_detect_dipoles(self):
        """Run detection on a difference image containing dipoles.
        """
        noiseLevel = 1.
        staticSeed = 1
        fluxLevel = 1000
        fluxRange = 1.5
        nSources = 10
        offset = 1
        dipoleFlag = "ip_diffim_DipoleFit_flag_classification"
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel, "fluxRange": fluxRange,
                  "nSrc": nSources}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        difference = science.clone()
        matchedTemplate.image.array[...] = np.roll(matchedTemplate.image.array[...], offset, axis=0)
        matchedTemplate.variance.array[...] = np.roll(matchedTemplate.variance.array[...], offset, axis=0)
        matchedTemplate.mask.array[...] = np.roll(matchedTemplate.mask.array[...], offset, axis=0)
        difference.maskedImage -= matchedTemplate.maskedImage[science.getBBox()]
        config = detectAndMeasure.DetectAndMeasureTask.ConfigClass()
        config.doApCorr = False
        config.doMerge = False
        config.doSkySources = False
        config.doForcedMeasurement = False
        detectionTask = detectAndMeasure.DetectAndMeasureTask(config=config)
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

        config.doMerge = True
        detectionTask2 = detectAndMeasure.DetectAndMeasureTask(config=config)
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
        noiseLevel = 1.
        staticSeed = 1
        transientSeed = 6
        transientFluxLevel = 1000.
        transientFluxRange = 1.5
        fluxLevel = 500
        kwargs = {"seed": staticSeed, "psfSize": 2.4, "fluxLevel": fluxLevel}
        science, sources = makeTestImage(noiseLevel=noiseLevel, noiseSeed=6, **kwargs)
        matchedTemplate, _ = makeTestImage(noiseLevel=noiseLevel/4, noiseSeed=7, **kwargs)
        config = detectAndMeasure.DetectAndMeasureTask.ConfigClass()
        config.doApCorr = False
        config.doMerge = False
        config.doSkySources = True
        config.doForcedMeasurement = False
        config.skySources.nSources = 5
        kernelWidth = np.max(science.psf.getKernel().getDimensions())//2
        detectionTask = detectAndMeasure.DetectAndMeasureTask(config=config)
        transients, transientSources = makeTestImage(seed=transientSeed, psfSize=2.4,
                                                     nSrc=10, fluxLevel=transientFluxLevel,
                                                     fluxRange=transientFluxRange,
                                                     noiseLevel=noiseLevel, noiseSeed=8)
        difference = science.clone()
        difference.maskedImage -= matchedTemplate.maskedImage
        difference.maskedImage += transients.maskedImage
        output = detectionTask.run(science, matchedTemplate, difference)
        skySources = output.diaSources[output.diaSources["sky_source"]]
        self.assertEqual(len(skySources), config.skySources.nSources)
        for skySource in skySources:
            # The sky sources should not be close to any other source
            with self.assertRaises(AssertionError):
                self._check_diaSource(transientSources, skySource, matchDistance=kernelWidth)
            with self.assertRaises(AssertionError):
                self._check_diaSource(sources, skySource, matchDistance=kernelWidth)
            # The sky sources should have low flux levels.
            self._check_diaSource(transientSources, skySource, matchDistance=1000, scale=0.,
                                  atol=np.sqrt(transientFluxRange*transientFluxLevel))

    def _check_diaSource(self, refSources, diaSource, refIds=None,
                         matchDistance=1., scale=1., usePsfFlux=True,
                         rtol=0.02, atol=None):
        """Match a diaSource with a source in a reference catalog
        and compare properties.
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


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
