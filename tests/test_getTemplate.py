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

import collections
import itertools
import unittest
import unittest.mock

import astropy.units as u
from astropy.stats import gaussian_sigma_to_fwhm
import numpy as np

import lsst.afw.cameraGeom.testUtils
import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.math
import lsst.afw.table
from lsst.daf.butler import DataCoordinate, DatasetRef, DatasetType, DimensionUniverse
import lsst.geom
import lsst.images
import lsst.images.psfs
from lsst.images.tests import compare_masked_image_to_legacy
import lsst.ip.diffim
import lsst.meas.algorithms
import lsst.meas.base.tests
import lsst.pex.exceptions
import lsst.pipe.base as pipeBase
import lsst.pipe.base.testUtils
import lsst.skymap
import lsst.utils.tests

from utils import generate_data_id, makeTestExposureRecord

# Change this to True, `setup display_ds9`, and open ds9 (or use another afw
# display backend) to show the tract/patch layouts on the image.
debug = False
if debug:
    import lsst.afw.display
    display = lsst.afw.display.Display()
    display.frame = 1

# Sum of squared warping kernel values along one axis, with the kernel
# normalized to unit sum at each sub-pixel phase, averaged over phase.
KERNEL_SUM_SQ_1D = {
    "bilinear": 2/3,
    "lanczos3": 0.8933699084,
    "lanczos4": 0.9141249535,
    "lanczos5": 0.9342245513,
}


def _showTemplate(box, template):
    """Show the corners of the template we made in this test."""
    for point in box.getCorners():
        display.dot("+", point.x, point.y, ctype="orange", size=40)
    display.frame = 2
    display.image(template, "warped template")
    display.frame = 3
    display.image(template.variance, "warped variance")


class GetTemplateTaskTestCase(lsst.utils.tests.TestCase):
    """Test that GetTemplateTask works on both one tract and multiple tract
    input coadd exposures.

    Makes a synthetic exposure large enough to fit four small tracts with 2x2
    (300x300 pixel) patches each, extracts pixels for those patches by warping,
    and tests GetTemplateTask's output against boxes that overlap various
    combinations of one or multiple tracts.
    """
    scale = 0.2  # arcsec/pixel
    # Use 5th order to minimize artifacts on the templates.
    coaddWarpKernel = "lanczos5"

    @classmethod
    def setUpClass(cls):
        # Warping the patches is the slowest part of every test, so it is done
        # once. Tests must not modify the exposure; the patch handles return
        # copies, and setUp gives each test its own lists of them.
        cls.skymap = cls._makeSkymap()
        cls.patches = collections.defaultdict(list)
        cls.dataIds = collections.defaultdict(list)
        cls.exposure = cls._makeExposure()

        if debug:
            display.image(cls.exposure, "base exposure")

        for tract_id in range(4):
            tract = cls.skymap.generateTract(tract_id)
            cls._makePatches(tract)

    @classmethod
    def tearDownClass(cls):
        del cls.skymap, cls.patches, cls.dataIds, cls.exposure

    def setUp(self):
        cls = type(self)
        self.patches = collections.defaultdict(list, {k: list(v) for k, v in cls.patches.items()})
        self.dataIds = collections.defaultdict(list, {k: list(v) for k, v in cls.dataIds.items()})
        self.coaddRefs = {}
        self.varianceBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))

    @classmethod
    def _makeSkymap(cls):
        """Make a Skymap with 4 tracts with 4 patches each.
        """
        tractScale = 0.02  # degrees
        # On-sky coordinates of the tract centers.
        coords = [(0, 0),
                  (0, tractScale),
                  (tractScale, 0),
                  (tractScale, tractScale),
                  ]
        config = lsst.skymap.DiscreteSkyMap.ConfigClass()
        config.raList = [c[0] for c in coords]
        config.decList = [c[1] for c in coords]
        # Half the tract center step size, to keep the tract overlap small.
        config.radiusList = [tractScale/2 for c in coords]
        config.projection = "TAN"
        config.pixelScale = cls.scale
        config.tractOverlap = 0.0005
        config.tractBuilder = "legacy"
        config.tractBuilder["legacy"].patchInnerDimensions = (300, 300)
        config.tractBuilder["legacy"].patchBorder = 10
        return lsst.skymap.DiscreteSkyMap(config=config)

    @classmethod
    def _makeExposure(cls):
        """Create a large image to break up into tracts and patches.

        The image will have a source every 100 pixels in x and y, and a WCS
        that results in the tracts all fitting in the image, with tract=0
        in the lower left, tract=1 to the right, tract=2 above, and tract=3
        to the upper right.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(-200, -200), lsst.geom.Point2I(800, 800))
        # This WCS was constructed so that tract 0 mostly fills the lower left
        # quadrant of the image, and the other tracts fill the rest; slight
        # extra rotation as a check on the final warp layout, scaled by 5%
        # from the patch pixel scale.
        cd_matrix = lsst.afw.geom.makeCdMatrix(1.05*cls.scale*lsst.geom.arcseconds, 93*lsst.geom.degrees)
        wcs = lsst.afw.geom.makeSkyWcs(lsst.geom.Point2D(120, 150),
                                       lsst.geom.SpherePoint(0, 0, lsst.geom.radians),
                                       cd_matrix)
        dataset = lsst.meas.base.tests.TestDataset(box, wcs=wcs)
        for x, y in itertools.product(np.arange(0, 500, 100), np.arange(0, 500, 100)):
            dataset.addSource(1e5, lsst.geom.Point2D(x, y))
        exposure, _ = dataset.realize(2, dataset.makeMinimalSchema())
        exposure.setFilter(lsst.afw.image.FilterLabel("a", "a_test"))
        return exposure

    @classmethod
    def _makePatches(cls, tract):
        """Populate the patches and dataId dicts, keyed on tract id, with the
        warps of the main exposure and minimal dataIds, respectively.
        """
        if debug:
            color = ['red', 'green', 'cyan', 'yellow'][tract.tract_id]
            point = cls.exposure.wcs.skyToPixel(tract.ctr_coord)
            # Show the tract center, colored by tract id.
            display.dot("x", point.x, point.y, ctype=color, size=30)

        config = lsst.afw.math.Warper.ConfigClass()
        config.warpingKernelName = cls.coaddWarpKernel
        warper = lsst.afw.math.Warper.fromConfig(config)
        for patchId in range(tract.num_patches.x*tract.num_patches.y):
            patch = tract.getPatchInfo(patchId)
            box = patch.getOuterBBox()

            if debug:
                # Show the patch corners as patch ids, colored by tract id.
                points = cls.exposure.wcs.skyToPixel(patch.wcs.pixelToSky([lsst.geom.Point2D(x)
                                                                           for x in box.getCorners()]))
                for p in points:
                    display.dot(patchId, p.x, p.y, ctype=color)

            # This is mostly taken from drp_tasks makePsfMatchedWarp, but
            # ip_diffim cannot depend on drp_tasks.
            xyTransform = lsst.afw.geom.makeWcsPairTransform(cls.exposure.wcs, patch.wcs)
            warpedPsf = lsst.meas.algorithms.WarpedPsf(cls.exposure.psf, xyTransform)
            warped = warper.warpExposure(patch.wcs, cls.exposure, destBBox=box)
            warped.setPsf(warpedPsf)

            warped.getInfo().setCoaddInputs(
                cls._makeCoaddInputs([(cls.exposure.wcs, cls.exposure.getBBox())]))
            dataRef = pipeBase.InMemoryDatasetHandle(
                warped,
                storageClass="ExposureF",
                copy=True,
                dataId=generate_data_id(
                    tract=tract,
                    patch=patch,
                )
            )
            cls.patches[tract.tract_id].append(dataRef)
            dataCoordinate = DataCoordinate.standardize({"tract": tract.tract_id,
                                                         "patch": patchId,
                                                         "band": "a",
                                                         "skymap": "skymap"},
                                                        universe=DimensionUniverse())
            cls.dataIds[tract.tract_id].append(dataCoordinate)

    def _checkMetadata(self, template, config, box, wcs, nPsfs):
        """Check that the various metadata components were set correctly.
        """
        expectedBox = lsst.geom.Box2I(box)
        expectedBox.grow(config.templateBorderSize)
        self.assertEqual(template.getBBox(), expectedBox)
        # WCS should match our exposure, not any of the coadd tracts.
        for tract in self.patches:
            self.assertNotEqual(template.wcs, self.patches[tract][0].get().wcs)
        self.assertEqual(template.wcs, self.exposure.wcs)
        # The template pixels are calibrated to nJy.
        self.assertEqual(template.photoCalib, lsst.afw.image.PhotoCalib(1.0))
        self.assertEqual(template.getXY0(), expectedBox.getMin())
        self.assertEqual(template.filter.bandLabel, "a")
        self.assertEqual(template.filter.physicalLabel, "a_test")
        # The template PSF is a Gaussian with the width of the CoaddPsf of
        # its inputs.
        self.assertIsInstance(template.psf, lsst.afw.detection.GaussianPsf)
        coaddPsf = self._makeCoaddPsf(template, config)
        self.assertEqual(coaddPsf.getComponentCount(), nPsfs)
        position = coaddPsf.getAveragePosition()
        self.assertFloatsAlmostEqual(template.psf.getSigma(),
                                     coaddPsf.computeShape(position).getTraceRadius())
        self.assertTrue(template.getInfo().hasCoaddInputs())
        self.assertEqual(len(template.getInfo().getCoaddInputs().ccds), nPsfs)

    def _makeConfig(self):
        """Make a GetTemplateTask config whose coaddWarpKernel is the kernel
        the test patches were made with.
        """
        config = lsst.ip.diffim.GetTemplateTask.ConfigClass()
        config.coaddWarpKernel = self.coaddWarpKernel
        return config

    def _makeCoaddPsf(self, template, config):
        """Return the CoaddPsf that ``run`` approximated with a Gaussian,
        rebuilt from the coadd inputs recorded on the template.
        """
        task = lsst.ip.diffim.GetTemplateTask(config=config)
        return task._makePsf(template, template.getInfo().getCoaddInputs().ccds, template.wcs)

    def _checkPixels(self, template, config, box):
        """Check that the pixel values in the template are close to the
        original image, calibrated to nJy.
        """
        # All pixels should have real values!
        expectedBox = lsst.geom.Box2I(box)
        expectedBox.grow(config.templateBorderSize)
        calibration = self.exposure.photoCalib.getCalibrationMean()
        expected = self.exposure.photoCalib.calibrateImage(self.exposure[expectedBox].maskedImage)

        if debug:
            _showTemplate(expectedBox, template)

        # Check that we fully filled the template from the patches.
        self.assertTrue(np.all(np.isfinite(template.image.array)))
        # Because of the scale changes, there will be some ringing in the
        # difference between the template and the original image; pick
        # tolerances large enough to account for that.
        self.assertImagesAlmostEqual(template.image, expected.image, rtol=.1, atol=4*calibration)
        # Variance plane ==4 in the original image (realize() takes a noise
        # sigma). Warping sets the level from the pixel areas and the two
        # warping kernels, which `_correctVariance` corrects. A per-pixel
        # ripple that no scalar correction can remove remains on top of that,
        # so check the level tightly and allow for the ripple around it.
        variance = template.variance.array[np.isfinite(template.variance.array)]
        median = np.median(variance)
        self.assertFloatsAlmostEqual(median, np.median(expected.variance.array),
                                     rtol=0.02, msg="variance level differs")
        self.assertLess(np.percentile(variance, 99)/median, 1.6, msg="variance ripple too large")
        self.assertGreater(np.percentile(variance, 1)/median, 0.6, msg="variance ripple too large")
        # Not checking the mask, as warping changes the sizes of the masks.

    def testRunOneTractInput(self):
        """Test a bounding box that fully fits inside one tract, with only
        that tract passed as input. This checks that the code handles a single
        tract input correctly.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        # Restrict to tract 0, since the box fits in just that tract.
        result = task.run(coaddExposureHandles={0: self.patches[0]},
                          bbox=box,
                          wcs=self.exposure.wcs,
                          dataIds={0: self.dataIds[0]},
                          physical_filter="a_test")

        # All 4 patches from tract 0 are included in this template.
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 4)
        self._checkPixels(result.template, task.config, box)

    def testRunOneTractMultipleInputs(self):
        """Test a bounding box that fully fits inside one tract but where
        multiple tracts were passed in. This checks that patches that are
        mostly NaN after warping are merged correctly in the output.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=box,
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")

        # All 4 patches from two tracts are included in this template.
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 6)
        self._checkPixels(result.template, task.config, box)

    def testRunTwoTracts(self):
        """Test a bounding box that crosses tract boundaries.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(200, 200), lsst.geom.Point2I(600, 600))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=box,
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")

        # All 4 patches from all 4 tracts are included in this template
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 9)
        self._checkPixels(result.template, task.config, box)

    def testRunCalibratesToNanojansky(self):
        """Test that the template is in nJy whatever the calibration of the
        coadds, and that coadds already in nJy are not rescaled.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        calibration = self.exposure.photoCalib.getCalibrationMean()
        self.assertNotEqual(calibration, 1.0)
        # The same coadd pixels, labelled as already in nJy.
        uncalibrated = []
        for handle in self.patches[0]:
            coadd = handle.get()
            coadd.setPhotoCalib(lsst.afw.image.PhotoCalib(1.0))
            uncalibrated.append(pipeBase.InMemoryDatasetHandle(coadd, storageClass="ExposureF", copy=True,
                                                               dataId=handle.dataId))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        templates = {}
        for key, handles in ((calibration, self.patches[0]), (1.0, uncalibrated)):
            result = task.run(coaddExposureHandles={0: handles},
                              bbox=box,
                              wcs=self.exposure.wcs,
                              dataIds={0: self.dataIds[0]},
                              physical_filter="a_test")
            self.assertEqual(result.template.photoCalib, lsst.afw.image.PhotoCalib(1.0))
            templates[key] = result.template

        self.assertFloatsAlmostEqual(templates[calibration].image.array,
                                     calibration*templates[1.0].image.array, rtol=1e-6)
        self.assertFloatsAlmostEqual(templates[calibration].variance.array,
                                     calibration**2*templates[1.0].variance.array, rtol=1e-6)

    def testGaussianPsfGridFallback(self):
        """If the CoaddPsf cannot be evaluated at its average position, the
        Gaussian width is averaged over a grid of positions instead.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        result = task.run(coaddExposureHandles={0: self.patches[0]},
                          bbox=box,
                          wcs=self.exposure.wcs,
                          dataIds={0: self.dataIds[0]},
                          physical_filter="a_test")
        template = result.template
        template.setPsf(self._makeCoaddPsf(template, task.config))
        expected = lsst.ip.diffim.utils.evaluateMeanPsfFwhm(
            template,
            fwhmExposureBuffer=task.config.fwhmExposureBuffer,
            fwhmExposureGrid=task.config.fwhmExposureGrid,
        )/gaussian_sigma_to_fwhm

        error = lsst.pex.exceptions.InvalidParameterError("No inputs at the average position.")
        with unittest.mock.patch("lsst.ip.diffim.getTemplate.getPsfFwhm", side_effect=error):
            with self.assertLogs(task.log.name, level="INFO") as cm:
                psf = task._makeGaussianPsf(template)

        self.assertIn("grid of points", "\n".join(cm.output))
        self.assertIsInstance(psf, lsst.afw.detection.GaussianPsf)
        self.assertFloatsAlmostEqual(psf.getSigma(), expected)
        self.assertFloatsAlmostEqual(task.metadata["templatePsfSigma"], expected)

    def testRunNoTemplate(self):
        """A bounding box that doesn't overlap the patches will raise.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(1200, 1200), lsst.geom.Point2I(1600, 1600))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        with self.assertRaisesRegex(lsst.pipe.base.NoWorkFound, "No patches found"):
            task.run(coaddExposureHandles=self.patches,
                     bbox=box,
                     wcs=self.exposure.wcs,
                     dataIds=self.dataIds,
                     physical_filter="a_test")

    def testMissingPatches(self):
        """Test that a missing patch results in an appropriate mask.

        This fixes the bug reported on DM-44997 (image and variance were NaN
        but the mask was not set to NO_DATA for those pixels).
        """
        # tract=0, patch=1 is the lower-left corner, as displayed in DS9.
        self.patches[0].pop(1)
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=box,
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        no_data = (result.template.mask.array & result.template.mask.getPlaneBitMask("NO_DATA")) != 0
        self.assertTrue(np.isfinite(result.template.image.array).all())
        self.assertTrue(np.isfinite(result.template.variance.array).all())
        self.assertEqual(no_data.sum(), 20990)

    @lsst.utils.tests.methodParameters(
        box=[
            lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180)),
            lsst.geom.Box2I(lsst.geom.Point2I(200, 200), lsst.geom.Point2I(600, 600)),
        ],
        nInput=[8, 16],
    )
    def testNanInputs(self, box=None, nInput=None):
        """Test that the template has finite values when some of the input
        pixels have NaN as variance.
        """
        for tract, patchRefs in self.patches.items():
            for patchRef in patchRefs:
                patchCoadd = patchRef.get()
                bbox = lsst.geom.Box2I()
                bbox.include(lsst.geom.Point2I(patchCoadd.getBBox().getCenter()))
                bbox.grow(3)
                patchCoadd.variance[bbox].array *= np.nan

        box = lsst.geom.Box2I(lsst.geom.Point2I(200, 200), lsst.geom.Point2I(600, 600))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=box,
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        if debug:
            _showTemplate(box, result.template)
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 9)
        # We just check that the pixel values are all finite. We cannot check that pixel values
        # in the template are closer to the original anymore.
        self.assertTrue(np.isfinite(result.template.image.array).all())

    def _runCorrection(self, doCorrectVarianceWarpCorrelation=False, doScaleVariance=False,
                       handles=None, instrument="SciCam", coaddWarpKernel=None,
                       templateWarpKernel="lanczos3"):
        """Run the task on tract 0, with a science image from ``instrument``.

        The plate scale correction is applied only when ``instrument`` is not
        the instrument of the template coadds, which is ``"TestCam"`` unless
        ``handles`` says otherwise. ``coaddWarpKernel`` defaults to the kernel
        the patches were made with, so the reconstructed factors are exact.
        """
        config = self._makeConfig()
        config.doCorrectVarianceWarpCorrelation = doCorrectVarianceWarpCorrelation
        config.doScaleVariance = doScaleVariance
        if coaddWarpKernel is not None:
            config.coaddWarpKernel = coaddWarpKernel
        config.warp.warpingKernelName = templateWarpKernel
        task = lsst.ip.diffim.GetTemplateTask(config=config)
        result = task.run(coaddExposureHandles={0: handles or self.patches[0]},
                          bbox=lsst.geom.Box2I(self.varianceBox),
                          wcs=self.exposure.wcs,
                          dataIds={0: self.dataIds[0]},
                          physical_filter="a_test",
                          instrument=instrument)
        return task, result.template

    def _makeScaledWcs(self, factor):
        """Make a WCS like the base exposure's, but with its pixel scale
        multiplied by ``factor``.
        """
        cdMatrix = lsst.afw.geom.makeCdMatrix(factor*1.05*self.scale*lsst.geom.arcseconds,
                                              93*lsst.geom.degrees)
        return lsst.afw.geom.makeSkyWcs(lsst.geom.Point2D(120, 150),
                                        lsst.geom.SpherePoint(0, 0, lsst.geom.radians),
                                        cdMatrix)

    @staticmethod
    def _makeCoaddInputs(records, instrument="TestCam"):
        """Make a CoaddInputs holding the given input records.

        Parameters
        ----------
        records : `list` [`tuple` [`lsst.afw.geom.SkyWcs` or `None`, \
                                   `lsst.geom.Box2I`]]
            The WCS and bbox to record for each input. A `None` WCS makes a
            record that `_plateScaleFactor` has to skip. May be empty, to
            simulate a coadd whose plate scale cannot be reconstructed.
        instrument : `str` or `None`, optional
            Instrument label to record for every input. `None` leaves the
            inputs without a visit info, and an empty string gives them a
            visit info with no label; either way the instrument is unknown.
        """
        ccdSchema = lsst.afw.table.ExposureTable.makeMinimalSchema()
        weightKey = ccdSchema.addField("weight", type=float, doc="Coadd weight")
        coaddInputs = lsst.afw.image.CoaddInputs(
            lsst.afw.table.ExposureTable.makeMinimalSchema(), ccdSchema)
        for wcs, bbox in records:
            record = coaddInputs.ccds.addNew()
            record.setWcs(wcs)
            record.setBBox(bbox)
            if instrument is not None:
                record.setVisitInfo(lsst.afw.image.VisitInfo(instrumentLabel=instrument))
            # Included because real coadds have it, though a single-record
            # correction does not use it.
            record.set(weightKey, 1.0)
        return coaddInputs

    def _patchHandles(self, tract, records, instrument="TestCam"):
        """Return handles for a tract's patches, with their CoaddInputs
        replaced by ``records`` from ``instrument``.
        """
        handles = []
        for ref in self.patches[tract]:
            coadd = ref.get()
            coadd.getInfo().setCoaddInputs(self._makeCoaddInputs(records, instrument=instrument))
            handles.append(pipeBase.InMemoryDatasetHandle(
                coadd, storageClass="ExposureF", copy=True, dataId=ref.dataId))
        return handles

    def testCorrectVarianceWarpCorrelation(self):
        """The correlation correction covers both warps: the one that built
        the coadds and the one this task performs.
        """
        _, off = self._runCorrection(instrument="TestCam")
        for coaddKernel, templateKernel in (("lanczos3", "lanczos3"),
                                            ("lanczos5", "lanczos3"),
                                            ("lanczos5", "lanczos5"),
                                            ("bilinear", "lanczos3")):
            with self.subTest(coaddWarpKernel=coaddKernel, templateWarpKernel=templateKernel):
                task, on = self._runCorrection(doCorrectVarianceWarpCorrelation=True,
                                               instrument="TestCam",
                                               coaddWarpKernel=coaddKernel,
                                               templateWarpKernel=templateKernel)
                factor = task.metadata["varianceWarpCorrelationFactor"]
                # A kernel spreads power, so this always makes the plane
                # larger.
                self.assertGreater(factor, 1.0)
                # The phase average converges slowly for bilinear, because
                # its sum of squares has a kink where the phase wraps.
                rtol = 1e-4 if "bilinear" in (coaddKernel, templateKernel) else 1e-5
                self.assertFloatsAlmostEqual(
                    factor,
                    1/(KERNEL_SUM_SQ_1D[coaddKernel]*KERNEL_SUM_SQ_1D[templateKernel])**2,
                    rtol=rtol)
                self.assertEqual(task.metadata["variancePlateScaleFactor"], 1.0)
                if templateKernel == "lanczos3":
                    # A pure scalar, with no per-pixel structure of its own;
                    # only comparable when `off` used the same template warp.
                    self.assertFloatsAlmostEqual(on.variance.array, off.variance.array*factor,
                                                 rtol=1e-5, ignoreNaNs=True)

    def testCorrectVariancePlateScale(self):
        """The plate scale correction is the total pixel area change from the
        images the coadds were built from to the science image.
        """
        # The science image and the coadds are from the same instrument, so
        # the factor is 1 without reading any coadd inputs.
        task, off = self._runCorrection(instrument="TestCam")
        self.assertEqual(task.metadata["variancePlateScaleFactor"], 1.0)
        self.assertNotIn("templateVarianceCorrectionFactor", task.metadata)

        # The fixture's records are the science image itself, so there is no
        # net change in pixel area even when the instruments differ.
        task, on = self._runCorrection()
        self.assertFloatsAlmostEqual(task.metadata["variancePlateScaleFactor"], 1.0, rtol=1e-6)
        self.assertNotIn("varianceWarpCorrelationFactor", task.metadata)

        # Coarser original pixels than science pixels, the DECam-template
        # case: the correction is the ratio of their areas.
        for scaleFactor in (1.315, 0.5):
            with self.subTest(scaleFactor=scaleFactor):
                handles = self._patchHandles(
                    0, [(self._makeScaledWcs(scaleFactor), self.exposure.getBBox())])
                task, on = self._runCorrection(handles=handles)
                factor = task.metadata["variancePlateScaleFactor"]
                self.assertFloatsAlmostEqual(factor, scaleFactor**2, rtol=1e-6)
                self.assertFloatsAlmostEqual(on.variance.array, off.variance.array*factor,
                                             rtol=1e-5, ignoreNaNs=True)

    def testVarianceCorrectionsCompose(self):
        """The two corrections are independent and multiply together.

        The patches are warped from the exposure and the template is warped
        back to its pixel grid, so the two warps cancel geometrically, leaving
        only the two kernels' correlation to undo. The substituted coadd
        inputs claim coarser original pixels than the exposure really has, so
        with both corrections applied the variance must come back to its
        original value times the claimed change in pixel area.
        """
        plateScale = 1.315
        handles = self._patchHandles(0, [(self._makeScaledWcs(plateScale), self.exposure.getBBox())])
        # A science image from the coadds' instrument applies no correction
        # at all; one from another instrument applies both.
        _, off = self._runCorrection(handles=handles, instrument="TestCam")
        task, on = self._runCorrection(doCorrectVarianceWarpCorrelation=True, handles=handles)

        # Both together is exactly the product of the two factors.
        self.assertFloatsAlmostEqual(task.metadata["variancePlateScaleFactor"], plateScale**2, rtol=1e-6)
        scale = (task.metadata["varianceWarpCorrelationFactor"]
                 * task.metadata["variancePlateScaleFactor"])
        self.assertFloatsAlmostEqual(task.metadata["templateVarianceCorrectionFactor"],
                                     scale, rtol=1e-12)
        self.assertFloatsAlmostEqual(on.variance.array, off.variance.array*scale,
                                     rtol=1e-5, ignoreNaNs=True)

        # Only the variance plane is touched. Nothing about the correction
        # feeds back into how the patches are combined.
        self.assertImagesEqual(on.image, off.image)
        self.assertMasksEqual(on.mask, off.mask)

        # realize() takes a noise sigma, so the original variance is 4, not 2,
        # before calibrating to nJy.
        calibration = self.exposure.photoCalib.getCalibrationMean()
        original = calibration**2*np.median(self.exposure.variance.array)
        self.assertFloatsAlmostEqual(np.median(on.variance.array), original*plateScale**2, rtol=0.02)

        # Which is a real correction: uncorrected, the plane is too small by
        # the two warping kernels' phase-averaged sums of squares.
        self.assertFloatsAlmostEqual(
            np.median(off.variance.array)/original,
            (KERNEL_SUM_SQ_1D[self.coaddWarpKernel]*KERNEL_SUM_SQ_1D["lanczos3"])**2, rtol=0.02)

    def testCorrectVariancePlateScaleUsesOneRecord(self):
        """Only the first usable coadd input is read.

        The spread between records is just the local pixel scale, a few
        tenths of a percent across a real focal plane, so one stands for all
        of them.
        """
        handles = self._patchHandles(
            0, [(None, self.exposure.getBBox()),
                (self._makeScaledWcs(2.0), self.exposure.getBBox()),
                (self.exposure.wcs, self.exposure.getBBox())])
        task, _ = self._runCorrection(handles=handles)

        # The first record with a WCS: the one without is skipped, and the
        # last (which would give 1.0) is never reached.
        self.assertFloatsAlmostEqual(task.metadata["variancePlateScaleFactor"], 4.0, rtol=1e-6)

    def testCorrectVariancePlateScaleSkipsDroppedCoadd(self):
        """The plate scale comes from the first coadd that is used, not from
        one dropped for its instrument.
        """
        bbox = self.exposure.getBBox()
        handles = (self._patchHandles(0, [(self._makeScaledWcs(2.0), bbox)], instrument="OtherCam")[:1]
                   + self._patchHandles(0, [(self._makeScaledWcs(1.315), bbox)])[1:])
        task, _ = self._runCorrection(handles=handles)
        self.assertFloatsAlmostEqual(task.metadata["variancePlateScaleFactor"], 1.315**2, rtol=1e-6)

    def testCorrectVariancePlateScaleNextTract(self):
        """If the first coadd kept in a tract cannot give the plate scale, the
        first coadd kept in the next tract is used.
        """
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        handles = {0: self._patchHandles(0, []),
                   1: self._patchHandles(1, [(self._makeScaledWcs(1.315), self.exposure.getBBox())])}
        task.run(coaddExposureHandles=handles,
                 bbox=lsst.geom.Box2I(self.varianceBox),
                 wcs=self.exposure.wcs,
                 dataIds={0: self.dataIds[0], 1: self.dataIds[1]},
                 physical_filter="a_test",
                 instrument="SciCam")
        self.assertFloatsAlmostEqual(task.metadata["variancePlateScaleFactor"], 1.315**2, rtol=1e-6)

    def testCorrectVariancePlateScaleNeedsCoaddInputs(self):
        """Without usable coadd inputs the plate scale cannot be
        reconstructed, so the template cannot be built.
        """
        handles = self._patchHandles(0, [])
        with self.assertRaisesRegex(RuntimeError, "cannot be corrected for the change in pixel area"):
            self._runCorrection(handles=handles)

    def _runTract0(self, handles):
        """Build a template from the given tract 0 patch handles.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        return task.run(coaddExposureHandles={0: handles},
                        bbox=box,
                        wcs=self.exposure.wcs,
                        dataIds={0: self.dataIds[0]},
                        physical_filter="a_test")

    def testRunRecordsInstrument(self):
        """The instrument of each coadd is recorded with the template's coadd
        inputs.
        """
        result = self._runTract0(self.patches[0])
        ccds = result.template.getInfo().getCoaddInputs().ccds
        self.assertEqual(len(ccds), 4)
        self.assertEqual({record["instrument"] for record in ccds}, {"TestCam"})

    def testRunUnknownInstrument(self):
        """A coadd whose inputs record no instrument is dropped when the
        other coadds record one, and kept, with the template's instrument
        ``"UNKNOWN"``, when none of them records one.
        """
        records = [(self.exposure.wcs, self.exposure.getBBox())]
        # Inputs with no visit info, and with a visit info but no label.
        for unknown in (None, ""):
            with self.subTest(label=unknown):
                handles = (self._patchHandles(0, records, instrument=unknown)[:1]
                           + self._patchHandles(0, records)[1:])
                ccds, dropped = self._runLogged({0: handles})
                self.assertEqual(len(ccds), 3)
                self.assertEqual({record["instrument"] for record in ccds}, {"TestCam"})
                self.assertEqual(len(dropped), 1)
                self.assertIn("made from UNKNOWN", dropped[0])

                # No coadd records an instrument.
                result = self._runTract0(self._patchHandles(0, records, instrument=unknown))
                ccds = result.template.getInfo().getCoaddInputs().ccds
                self.assertEqual(len(ccds), 4)
                self.assertEqual({record["instrument"] for record in ccds}, {"UNKNOWN"})

    def _runLogged(self, handles, dataIds=None, box=None, instrument=None):
        """Build a template, capturing the messages about dropped coadds.

        Parameters
        ----------
        handles : `dict` [`int`, `list`]
            Patch handles, keyed by tract.
        dataIds : `dict` [`int`, `list`], optional
            Data ids matching ``handles``; defaults to the fixture's.
        box : `lsst.geom.Box2I`, optional
            Box to build the template on; defaults to one in tract 0.
        instrument : `str`, optional
            Instrument of the science image.

        Returns
        -------
        ccds : `lsst.afw.table.ExposureCatalog`
            The coadd inputs recorded on the template.
        dropped : `list` [`str`]
            The info messages about dropped coadds.
        """
        if dataIds is None:
            dataIds = {tract: self.dataIds[tract] for tract in handles}
        if box is None:
            box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask(config=self._makeConfig())
        with self.assertLogs(task.log.name, level="INFO") as cm:
            result = task.run(coaddExposureHandles=handles,
                              bbox=box,
                              wcs=self.exposure.wcs,
                              dataIds=dataIds,
                              physical_filter="a_test",
                              instrument=instrument)
        dropped = [record for record in cm.records
                   if record.getMessage().startswith("Dropping template coadd")]
        self.assertEqual({record.levelname for record in dropped}, {"INFO"} if dropped else set())
        return result.template.getInfo().getCoaddInputs().ccds, [record.getMessage() for record in dropped]

    def testRunDropsMinorityInstrument(self):
        """Coadds made from an instrument other than the one most coadds were
        made from are dropped.
        """
        records = [(self.exposure.wcs, self.exposure.getBBox())]
        handles = (self._patchHandles(0, records)[:3]
                   + self._patchHandles(0, records, instrument="OtherCam")[3:])

        ccds, dropped = self._runLogged({0: handles})
        self.assertEqual({record["patch"] for record in ccds}, {0, 1, 2})
        self.assertEqual({record["instrument"] for record in ccds}, {"TestCam"})
        self.assertEqual(len(dropped), 1)
        self.assertIn(str(self.dataIds[0][3]), dropped[0])
        self.assertIn("made from OtherCam", dropped[0])

    def testRunPrefersScienceInstrument(self):
        """Coadds made from the science image's instrument are used even when
        they are in the minority; otherwise the majority instrument is used.
        """
        records = [(self.exposure.wcs, self.exposure.getBBox())]
        handles = (self._patchHandles(0, records, instrument="OtherCam")[:3]
                   + self._patchHandles(0, records)[3:])

        ccds, dropped = self._runLogged({0: handles}, instrument="TestCam")
        self.assertEqual({record["patch"] for record in ccds}, {3})
        self.assertEqual({record["instrument"] for record in ccds}, {"TestCam"})
        self.assertEqual(len(dropped), 3)
        self.assertTrue(all("made from OtherCam" in message for message in dropped))

        # No coadd was made from the science image's instrument.
        ccds, dropped = self._runLogged({0: handles}, instrument="SciCam")
        self.assertEqual({record["patch"] for record in ccds}, {0, 1, 2})
        self.assertEqual({record["instrument"] for record in ccds}, {"OtherCam"})
        self.assertEqual(len(dropped), 1)

    def testRunInstrumentTie(self):
        """A tie goes to the instrument whose name sorts first.
        """
        records = [(self.exposure.wcs, self.exposure.getBBox())]
        handles = (self._patchHandles(0, records)[:2]
                   + self._patchHandles(0, records, instrument="OtherCam")[2:])

        ccds, dropped = self._runLogged({0: handles})
        self.assertEqual({record["patch"] for record in ccds}, {2, 3})
        self.assertEqual({record["instrument"] for record in ccds}, {"OtherCam"})
        self.assertEqual(len(dropped), 2)

    def testRunDropsNonOverlappingCoadd(self):
        """A coadd from another instrument is dropped even if it does not
        overlap the template.
        """
        records = [(self.exposure.wcs, self.exposure.getBBox())]
        # Tract 3 does not overlap a box in the lower left of tract 0.
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(50, 50))
        handles = {0: self.patches[0], 3: self._patchHandles(3, records, instrument="OtherCam")[:1]}
        dataIds = {0: self.dataIds[0], 3: self.dataIds[3][:1]}

        ccds, dropped = self._runLogged(handles, dataIds=dataIds, box=box)
        self.assertEqual({record["tract"] for record in ccds}, {0})
        self.assertEqual({record["instrument"] for record in ccds}, {"TestCam"})
        self.assertEqual(len(dropped), 1)
        self.assertIn(str(self.dataIds[3][0]), dropped[0])

    def _scaleInputVariance(self, tract, factor):
        """Return fresh handles for one tract's patches, with their variance
        planes multiplied by ``factor``.

        Parameters
        ----------
        tract : `int`
            Id of the tract whose patches should be copied.
        factor : `float`
            Factor to multiply the input variance planes by.

        Returns
        -------
        handles : `list` [`lsst.pipe.base.InMemoryDatasetHandle`]
            Handles to the modified patches.
        """
        handles = []
        for handle in self.patches[tract]:
            # ``copy=True`` on the original handles means this is a copy, so
            # the patches shared with the other tests are left untouched.
            patch = handle.get()
            patch.variance.array *= factor
            handles.append(pipeBase.InMemoryDatasetHandle(patch,
                                                          storageClass="ExposureF",
                                                          copy=True,
                                                          dataId=handle.dataId))
        return handles

    def testScaleVariance(self):
        """Test that the template variance plane is rescaled to match the
        empirical pixel noise, and that the factor used is recorded in the
        task metadata.
        """
        scaleFactor = 1.345
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))

        def _configureAndRunTask(doScaleVariance, varianceScale=1.):
            """Build a template from tract 0, optionally rescaling the input
            variance planes by ``varianceScale`` first.
            """
            config = self._makeConfig()
            config.doScaleVariance = doScaleVariance
            task = lsst.ip.diffim.GetTemplateTask(config=config)
            result = task.run(coaddExposureHandles={0: self._scaleInputVariance(0, varianceScale)},
                              bbox=box,
                              wcs=self.exposure.wcs,
                              dataIds={0: self.dataIds[0]},
                              physical_filter="a_test")
            return task, result.template

        # With scaling disabled the subtask is never constructed, and nothing
        # is recorded in the metadata.
        taskOff, templateOff = _configureAndRunTask(doScaleVariance=False)
        self.assertFalse(hasattr(taskOff, "scaleVariance"))
        self.assertNotIn("scaleTemplateVarianceFactor", taskOff.metadata)

        # The image is correlated after warping, so ``doScaleVariance`` will
        # apply a correction. The value here is empirical from a warped image
        # with white noise and no sources.
        taskOn, templateOn = _configureAndRunTask(doScaleVariance=True)
        factor = taskOn.metadata["scaleTemplateVarianceFactor"]
        self.assertFloatsAlmostEqual(factor, 1.145, atol=0.01,
                                     msg="Measured template variance scaling changed; check the"
                                         " warping kernels used by _makePatches and the task,"
                                         " and the order of the variance corrections.")
        # The only difference from the unscaled template is the constant
        # factor applied to the variance plane.
        self.assertFloatsAlmostEqual(templateOn.variance.array,
                                     templateOff.variance.array*factor, rtol=1e-5)
        # Tolerance here is float32 round-off: repeated runs of the task are
        # not bitwise identical.
        self.assertImagesAlmostEqual(templateOn.image, templateOff.image, rtol=1e-5, atol=1e-5)

        # If the input variance planes under-estimate the noise by a known
        # factor, the measured factor grows by that amount and the same
        # output variance plane is recovered.
        taskLow, templateLow = _configureAndRunTask(doScaleVariance=True, varianceScale=1/scaleFactor)
        self.assertFloatsAlmostEqual(taskLow.metadata["scaleTemplateVarianceFactor"],
                                     factor*scaleFactor, rtol=1e-5)
        self.assertImagesAlmostEqual(templateLow.variance, templateOn.variance, rtol=1e-5)

    def _runLegacyForFuture(self, raiseOnUndefinedMaskMap=True, handles=None):
        """Build a template from tract 0 with a task configured for the
        future output type, without converting it.

        The task is kept as ``self.futureTask``.

        Parameters
        ----------
        raiseOnUndefinedMaskMap : `bool`, optional
            Value of the task config field of the same name.
        handles : `list` [`lsst.pipe.base.InMemoryDatasetHandle`], optional
            Handles to the tract 0 patches to use instead of the fixture's.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            The legacy output struct.
        box : `lsst.geom.Box2I`
            The bounding box the template was requested on, before the task
            grew it by the template border.
        """
        config = self._makeConfig()
        config.output_image_type = "future"
        config.raiseOnUndefinedMaskMap = raiseOnUndefinedMaskMap
        task = lsst.ip.diffim.GetTemplateTask(config=config)
        self.visit = 9876
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        result = task.run(coaddExposureHandles={0: self.patches[0] if handles is None else handles},
                          bbox=box,
                          wcs=self.exposure.wcs,
                          dataIds={0: self.dataIds[0]},
                          physical_filter="a_test",
                          visit=self.visit)
        # run() always produces legacy types, whatever output_image_type is.
        self.assertIsInstance(result.template, lsst.afw.image.ExposureF)
        self.futureTask = task
        return result, box

    def _runFuture(self, handles=None):
        """Build a template and convert it, the way ``runQuantum`` does.

        Parameters
        ----------
        handles : `list` [`lsst.pipe.base.InMemoryDatasetHandle`], optional
            Handles to the tract 0 patches to use instead of the fixture's.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            The converted output struct.
        box : `lsst.geom.Box2I`
            The bounding box the template was requested on, before the task
            grew it by the template border.
        """
        # lsst.images needs per-amplifier raw geometry and a field angle
        # transform, which the trivial test detectors do not have.
        detector = list(lsst.afw.cameraGeom.testUtils.CameraWrapper().camera)[0]
        result, box = self._runLegacyForFuture(handles=handles)
        result.template = self.futureTask.convert_outputs_to_future(
            result.template, self._coaddRefs(0), detector=detector,
            exposureRecord=self._exposureRecord())
        return result, box, detector

    def _exposureRecord(self):
        """Return the exposure record that ``runQuantum`` rebuilds from the
        observation metadata of the science image.
        """
        return makeTestExposureRecord(DimensionUniverse(), visit=self.visit)

    def _coaddRefs(self, tract):
        """Return butler references for the coadds of one tract, like the
        ones ``runQuantum`` passes on from its ``coaddExposures`` input.

        Parameters
        ----------
        tract : `int`
            Id of the tract whose coadds to make references for.

        Returns
        -------
        refs : `list` [`lsst.daf.butler.DatasetRef`]
            One reference per patch of that tract. The same references are
            returned on every call, because each new one would get a new
            dataset id.
        """
        if tract not in self.coaddRefs:
            datasetType = DatasetType("template_coadd", ("tract", "patch", "band", "skymap"),
                                      "ExposureF", universe=DimensionUniverse())
            self.coaddRefs[tract] = [DatasetRef(datasetType, dataId, run="test_run")
                                     for dataId in self.dataIds[tract]]
        return self.coaddRefs[tract]

    def testConvertOutputsToFuture(self):
        """Test that the output template is converted to a DifferenceImage
        that keeps the pixels and the detector it was built on.
        """
        detector = list(lsst.afw.cameraGeom.testUtils.CameraWrapper().camera)[0]
        result, box = self._runLegacyForFuture()
        legacy = result.template.clone()
        result.template = self.futureTask.convert_outputs_to_future(
            result.template, self._coaddRefs(0), detector=detector,
            exposureRecord=self._exposureRecord())

        self.assertIsInstance(result.template, lsst.images.DifferenceImage)
        compare_masked_image_to_legacy(result.template, legacy.maskedImage,
                                       plane_map=lsst.images.get_legacy_template_mask_planes())
        self.assertEqual(result.template.unit, u.nJy)
        self.assertIsNotNone(result.template.detector)
        self.assertEqual(result.template.detector.name, detector.getName())
        # The template is grown by the border, so it is larger than both the
        # requested box and the detector.
        expectedBox = lsst.geom.Box2I(box)
        expectedBox.grow(lsst.ip.diffim.GetTemplateTask.ConfigClass().templateBorderSize)
        self.assertEqual(result.template.bbox.to_legacy(), expectedBox)
        self.assertEqual(result.template.obs_info.visit_id, self.visit)

    def testConvertOutputsToFutureCoaddMaskPlanes(self):
        """Test that the coadd mask planes a warped template carries are
        converted, and come back when the image is converted to legacy.

        The template is converted with the template plane map, so the coadd
        planes are named there; the source injection planes are optional, and
        are added because this template has pixels set in them.
        """
        detector = list(lsst.afw.cameraGeom.testUtils.CameraWrapper().camera)[0]
        planes = ("CLIPPED", "REJECTED", "INEXACT_PSF", "SENSOR_EDGE", "HIGH_VARIANCE",
                  "INJECTED", "INJECTED_CORE")
        result, _ = self._runLegacyForFuture()
        mask = result.template.mask
        for n, plane in enumerate(planes):
            mask.addMaskPlane(plane)
            mask.array[0, n] |= mask.getPlaneBitMask(plane)
        legacyTemplate = result.template.clone()

        result.template = self.futureTask.convert_outputs_to_future(
            result.template, self._coaddRefs(0), detector=detector,
            exposureRecord=self._exposureRecord())

        template = result.template
        self.assertIsInstance(template, lsst.images.DifferenceImage)
        # Name the optional planes that have pixels set in the map, so that
        # dropping one fails the comparison instead of skipping it. The others
        # are left out because they may be registered by unrelated tests.
        optional = lsst.images.get_legacy_optional_mask_planes()
        planeMap = lsst.images.get_legacy_template_mask_planes()
        planeMap.update({plane: optional[plane] for plane in planes if plane in optional})
        compare_masked_image_to_legacy(template, legacyTemplate.maskedImage, plane_map=planeMap)
        # The butler converts a DifferenceImage to an ExposureF with the
        # difference image plane map, which does not name the coadd planes;
        # they keep their own names instead of being dropped.
        legacy = template.to_legacy()
        for n, plane in enumerate(planes):
            with self.subTest(plane=plane):
                self.assertIn(plane, legacy.mask.getMaskPlaneDict())
                self.assertEqual(
                    np.count_nonzero(legacy.mask.array & legacy.mask.getPlaneBitMask(plane)), 1)

    def _addUnmappedMaskPlane(self, template):
        """Set one pixel in the template's BRIGHT_OBJECT plane, which the
        template plane map does not include, and one in its CLIPPED plane,
        which it does.
        """
        mask = template.mask
        for n, plane in enumerate(("BRIGHT_OBJECT", "CLIPPED")):
            mask.addMaskPlane(plane)
            mask.array[0, n] |= mask.getPlaneBitMask(plane)

    def testConvertOutputsToFutureUnmappedMaskPlaneRaises(self):
        """With raiseOnUndefinedMaskMap=True, a template with pixels set in
        an unmapped mask plane cannot be converted.
        """
        detector = list(lsst.afw.cameraGeom.testUtils.CameraWrapper().camera)[0]
        result, _ = self._runLegacyForFuture()
        self._addUnmappedMaskPlane(result.template)
        with self.assertRaisesRegex(RuntimeError, "BRIGHT_OBJECT"):
            self.futureTask.convert_outputs_to_future(
                result.template, self._coaddRefs(0), detector=detector,
                exposureRecord=self._exposureRecord())

    def testConvertOutputsToFutureUnmappedMaskPlaneWarns(self):
        """With raiseOnUndefinedMaskMap=False, an unmapped mask plane is
        dropped with a warning and the mapped planes are kept.
        """
        detector = list(lsst.afw.cameraGeom.testUtils.CameraWrapper().camera)[0]
        result, _ = self._runLegacyForFuture(raiseOnUndefinedMaskMap=False)
        self._addUnmappedMaskPlane(result.template)
        with self.assertLogs(self.futureTask.log.name, level="WARNING") as cm:
            template = self.futureTask.convert_outputs_to_future(
                result.template, self._coaddRefs(0), detector=detector,
                exposureRecord=self._exposureRecord())
        self.assertIn("BRIGHT_OBJECT", "\n".join(cm.output))
        self.assertIsInstance(template, lsst.images.DifferenceImage)
        self.assertNotIn("BRIGHT_OBJECT", template.mask.schema.names)
        self.assertEqual(np.count_nonzero(template.mask.get("CLIPPED")), 1)

    def testConvertOutputsToFutureLosesProvenance(self):
        """The coadd inputs that `run` attaches are not carried by
        `lsst.images.DifferenceImage`.
        """
        legacy, _ = self._runLegacyForFuture()
        self.assertTrue(legacy.template.getInfo().hasCoaddInputs())
        result, _, _ = self._runFuture()
        self.assertFalse(result.template.to_legacy().getInfo().hasCoaddInputs())

    def testConvertOutputsToFutureTemplates(self):
        """The coadds that went into the template are recorded on the
        converted image.
        """
        result, _, _ = self._runFuture()
        refs = {(ref.dataId["tract"], ref.dataId["patch"]): ref for ref in self._coaddRefs(0)}

        templates = result.template.templates
        self.assertGreater(len(templates), 0)
        for template in templates:
            with self.subTest(tract=template.tract, patch=template.patch):
                ref = refs[template.tract, template.patch]
                self.assertEqual(template.skymap, "skymap")
                self.assertEqual(template.dataset_id, ref.id)
                self.assertEqual(template.dataset_run, ref.run)
                self.assertEqual(template.instrument, "TestCam")
                self.assertFalse(template.psf_shape_flag)
                self.assertGreater(template.psf_shape_xx, 0)
                self.assertGreater(template.bounds.area, 0)

    def testConvertOutputsToFutureInstrument(self):
        """The template instrument is the one recorded by the coadds that
        know it, or ``"UNKNOWN"`` if none of them do.
        """
        records = [(self.exposure.wcs, self.exposure.getBBox())]
        handles = (self._patchHandles(0, records, instrument=None)[:1]
                   + self._patchHandles(0, records)[1:])
        result, _, _ = self._runFuture(handles=handles)
        self.assertEqual({template.instrument for template in result.template.templates}, {"TestCam"})

        result, _, _ = self._runFuture(handles=self._patchHandles(0, records, instrument=None))
        self.assertEqual({template.instrument for template in result.template.templates}, {"UNKNOWN"})

    def testConvertOutputsToFuturePsf(self):
        """The Gaussian PSF set by `run` is kept, defined over the detector.
        """
        result, _ = self._runLegacyForFuture()
        detector = list(lsst.afw.cameraGeom.testUtils.CameraWrapper().camera)[0]
        legacyPsf = result.template.getPsf()

        result.template = self.futureTask.convert_outputs_to_future(
            result.template, self._coaddRefs(0), detector=detector,
            exposureRecord=self._exposureRecord())

        psf = result.template.psf
        self.assertIsInstance(psf, lsst.images.psfs.GaussianPointSpreadFunction)
        self.assertFloatsAlmostEqual(psf.sigma, legacyPsf.getSigma())
        self.assertEqual(psf.bounds.to_legacy(), detector.getBBox())
        self.assertEqual(psf.kernel_bbox.to_legacy().getDimensions(), legacyPsf.getDimensions())

    def testConvertOutputsToFuturePsfToLegacy(self):
        """The Gaussian survives conversion back to an Exposure, which is how
        a task that has not been converted reads the template.
        """
        result, _, _ = self._runFuture()

        psf = result.template.to_legacy().getPsf()
        self.assertIsInstance(psf, lsst.afw.detection.GaussianPsf)
        self.assertFloatsAlmostEqual(psf.getSigma(), result.template.psf.sigma)


class WarpCorrelationFactorTestCase(lsst.utils.tests.TestCase):
    """Test the warping correlation correction on pure noise, with no other
    variance correction.
    """

    def setUp(self):
        self.size = 1024
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(self.size, self.size))

    def _makeWcs(self, angle, offset):
        """Make a WCS with a fixed plate scale, rotated by ``angle`` degrees
        and with its reference pixel moved by ``offset`` pixels.
        """
        center = lsst.geom.Point2D(self.size/2 + offset[0], self.size/2 + offset[1])
        cdMatrix = lsst.afw.geom.makeCdMatrix(0.2*lsst.geom.arcseconds, angle*lsst.geom.degrees)
        return lsst.afw.geom.makeSkyWcs(center, lsst.geom.SpherePoint(45, 30, lsst.geom.degrees),
                                        cdMatrix)

    def testCorrelationFactorRestoresNoisePower(self):
        """Warp pure noise to a rotated, offset frame with the coadd kernel,
        then to a third frame with the task's own warper.

        Warping leaves the variance of a flux summed over a region (the noise
        power) unchanged, but lowers the per-pixel scatter and the variance
        plane. The variance plane is lowered by a larger factor because the
        input to the second warp is already correlated. The correction must
        bring the variance plane back to the noise power of the original image.
        """
        sigma = 2.0
        rng = np.random.default_rng(42)
        original = lsst.afw.image.ExposureF(self.bbox, self._makeWcs(0, (0, 0)))
        original.image.array[:] = rng.normal(0, sigma, original.image.array.shape)
        original.variance.array[:] = sigma**2

        config = lsst.ip.diffim.GetTemplateTask.ConfigClass()
        config.doScaleVariance = False
        config.coaddWarpKernel = "lanczos3"
        task = lsst.ip.diffim.GetTemplateTask(config=config)
        coadd = lsst.afw.math.Warper(config.coaddWarpKernel).warpExposure(
            self._makeWcs(37, (0.31, 0.72)), original, destBBox=self.bbox)
        template = task.warper.warpExposure(self._makeWcs(-21, (0.58, 0.23)), coadd, destBBox=self.bbox)

        # The central quarter is covered in every frame.
        inner = lsst.geom.Box2I(lsst.geom.Point2I(self.size//4, self.size//4),
                                lsst.geom.Extent2I(self.size//2, self.size//2))
        originalStd = original.image[inner].array.std()
        image = template.image[inner].array.astype(float)
        pixelStd = image.std()
        # Each warp moves flux between neighboring pixels but doesn't create
        # or destroy any, so the noise in a flux summed over a region is the
        # same before and after. Summing the image over blocks and measuring
        # how much those sums scatter gives that quantity directly, without
        # depending on the variance plane. Every placement of the block is
        # used, from a summed-area table, to reduce the noise.
        blockSize = 8
        table = np.pad(image.cumsum(axis=0).cumsum(axis=1), ((1, 0), (1, 0)))
        sums = (table[blockSize:, blockSize:] - table[:-blockSize, blockSize:]
                - table[blockSize:, :-blockSize] + table[:-blockSize, :-blockSize])
        blockStd = sums.std()/blockSize
        planeStdBefore = np.sqrt(template.variance[inner].array.mean())
        task._correctVariance(template, plateScaleFactor=1.0)
        planeStdAfter = np.sqrt(template.variance[inner].array.mean())
        factor = task.metadata["varianceWarpCorrelationFactor"]

        # afw propagates the variance as if each warp's input pixels were
        # independent, so the plane drops by exactly the factor's reciprocal.
        self.assertFloatsAlmostEqual(planeStdBefore, sigma/np.sqrt(factor), rtol=1e-5)
        self.assertFloatsAlmostEqual(planeStdAfter, sigma, rtol=1e-5)

        # The warps conserve the noise power, which the corrected plane
        # matches.
        self.assertFloatsAlmostEqual(blockStd, originalStd, rtol=0.03)
        self.assertFloatsAlmostEqual(blockStd, planeStdAfter, rtol=0.03)

        # The per-pixel scatter drops, though by less than the uncorrected
        # plane, because the second warp receives correlated input. The
        # corrected plane deliberately exceeds it.
        self.assertLess(pixelStd, 0.9*originalStd)
        self.assertGreater(pixelStd, planeStdBefore)
        self.assertLess(pixelStd, planeStdAfter)


class GetTemplateConnectionsTestCase(lsst.utils.tests.TestCase):
    """Test the connections that ``output_image_type`` switches on.

    These only need the config classes, not a built template.
    """

    def testOutputImageType(self):
        """Test that the future output type changes the template storage
        class, and leaves the input connections alone.
        """
        Connections = lsst.ip.diffim.GetTemplateTask.ConfigClass.ConnectionsClass
        config = lsst.ip.diffim.GetTemplateTask.ConfigClass()

        connections = Connections(config=config)
        self.assertEqual(connections.template.storageClass, "ExposureF")
        self.assertNotIn("obs_info", connections.inputs)

        config.output_image_type = "future"
        connections = Connections(config=config)
        self.assertEqual(connections.template.storageClass, "DifferenceImage")
        # The dataset name is unchanged; only its storage class differs.
        self.assertEqual(connections.template.name, "goodSeeingDiff_templateExp")
        # The detector the conversion needs comes from the science image.
        self.assertEqual(connections.detector.name, "calexp.detector")
        self.assertEqual(connections.detector.storageClass, "Detector")
        # The observation metadata is only read in future mode, where it is
        # the one input the conversion adds.
        self.assertEqual(connections.obs_info.name, "calexp.obs_info")
        self.assertEqual(connections.obs_info.storageClass, "ObservationInfo")

    def testLintConnections(self):
        """Check that the connections are self-consistent in both modes.
        """
        for task in (lsst.ip.diffim.GetTemplateTask, lsst.ip.diffim.GetDcrTemplateTask):
            with self.subTest(task=task.__name__):
                lsst.pipe.base.testUtils.lintConnections(task.ConfigClass.ConnectionsClass)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
