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
    def setUp(self):
        self.scale = 0.2  # arcsec/pixel
        self.skymap = self._makeSkymap()
        self.patches = collections.defaultdict(list)
        self.dataIds = collections.defaultdict(list)
        self.coaddRefs = {}
        self.exposure = self._makeExposure()
        self.varianceBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))

        if debug:
            display.image(self.exposure, "base exposure")

        for tract_id in range(4):
            tract = self.skymap.generateTract(tract_id)
            self._makePatches(tract)

    def _makeSkymap(self):
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
        config.pixelScale = self.scale
        config.tractOverlap = 0.0005
        config.tractBuilder = "legacy"
        config.tractBuilder["legacy"].patchInnerDimensions = (300, 300)
        config.tractBuilder["legacy"].patchBorder = 10
        return lsst.skymap.DiscreteSkyMap(config=config)

    def _makeExposure(self):
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
        cd_matrix = lsst.afw.geom.makeCdMatrix(1.05*self.scale*lsst.geom.arcseconds, 93*lsst.geom.degrees)
        wcs = lsst.afw.geom.makeSkyWcs(lsst.geom.Point2D(120, 150),
                                       lsst.geom.SpherePoint(0, 0, lsst.geom.radians),
                                       cd_matrix)
        dataset = lsst.meas.base.tests.TestDataset(box, wcs=wcs)
        for x, y in itertools.product(np.arange(0, 500, 100), np.arange(0, 500, 100)):
            dataset.addSource(1e5, lsst.geom.Point2D(x, y))
        exposure, _ = dataset.realize(2, dataset.makeMinimalSchema())
        exposure.setFilter(lsst.afw.image.FilterLabel("a", "a_test"))
        return exposure

    def _makePatches(self, tract):
        """Populate the patches and dataId dicts, keyed on tract id, with the
        warps of the main exposure and minimal dataIds, respectively.
        """
        if debug:
            color = ['red', 'green', 'cyan', 'yellow'][tract.tract_id]
            point = self.exposure.wcs.skyToPixel(tract.ctr_coord)
            # Show the tract center, colored by tract id.
            display.dot("x", point.x, point.y, ctype=color, size=30)

        # Use 5th order to minimize artifacts on the templates.
        config = lsst.afw.math.Warper.ConfigClass()
        config.warpingKernelName = "lanczos5"
        warper = lsst.afw.math.Warper.fromConfig(config)
        for patchId in range(tract.num_patches.x*tract.num_patches.y):
            patch = tract.getPatchInfo(patchId)
            box = patch.getOuterBBox()

            if debug:
                # Show the patch corners as patch ids, colored by tract id.
                points = self.exposure.wcs.skyToPixel(patch.wcs.pixelToSky([lsst.geom.Point2D(x)
                                                                           for x in box.getCorners()]))
                for p in points:
                    display.dot(patchId, p.x, p.y, ctype=color)

            # This is mostly taken from drp_tasks makePsfMatchedWarp, but
            # ip_diffim cannot depend on drp_tasks.
            xyTransform = lsst.afw.geom.makeWcsPairTransform(self.exposure.wcs, patch.wcs)
            warpedPsf = lsst.meas.algorithms.WarpedPsf(self.exposure.psf, xyTransform)
            warped = warper.warpExposure(patch.wcs, self.exposure, destBBox=box)
            warped.setPsf(warpedPsf)

            warped.getInfo().setCoaddInputs(
                self._makeCoaddInputs([(self.exposure.wcs, self.exposure.getBBox())]))
            dataRef = pipeBase.InMemoryDatasetHandle(
                warped,
                storageClass="ExposureF",
                copy=True,
                dataId=generate_data_id(
                    tract=tract,
                    patch=patch,
                )
            )
            self.patches[tract.tract_id].append(dataRef)
            dataCoordinate = DataCoordinate.standardize({"tract": tract.tract_id,
                                                         "patch": patchId,
                                                         "band": "a",
                                                         "skymap": "skymap"},
                                                        universe=DimensionUniverse())
            self.dataIds[tract.tract_id].append(dataCoordinate)

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
        # warping kernels, which `_correctVariance` corrects up to the
        # difference between the configured coaddWarpKernel and the lanczos5
        # `_makePatches` really used. A per-pixel ripple that no scalar
        # correction can remove remains on top of that, so check the level
        # tightly and allow for the ripple around it.
        variance = template.variance.array[np.isfinite(template.variance.array)]
        median = np.median(variance)
        self.assertFloatsAlmostEqual(median, np.median(expected.variance.array),
                                     rtol=0.35, msg="variance level differs")
        self.assertLess(np.percentile(variance, 99)/median, 1.6, msg="variance ripple too large")
        self.assertGreater(np.percentile(variance, 1)/median, 0.6, msg="variance ripple too large")
        # Not checking the mask, as warping changes the sizes of the masks.

    def testRunOneTractInput(self):
        """Test a bounding box that fully fits inside one tract, with only
        that tract passed as input. This checks that the code handles a single
        tract input correctly.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask()
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
        task = lsst.ip.diffim.GetTemplateTask()
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
        task = lsst.ip.diffim.GetTemplateTask()
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
        task = lsst.ip.diffim.GetTemplateTask()
        templates = {}
        for photoCalib in (self.exposure.photoCalib, lsst.afw.image.PhotoCalib(1.0)):
            self._remakePatches(photoCalib)
            result = task.run(coaddExposureHandles={0: self.patches[0]},
                              bbox=box,
                              wcs=self.exposure.wcs,
                              dataIds={0: self.dataIds[0]},
                              physical_filter="a_test")
            self.assertEqual(result.template.photoCalib, lsst.afw.image.PhotoCalib(1.0))
            templates[photoCalib.getCalibrationMean()] = result.template

        self.assertFloatsAlmostEqual(templates[calibration].image.array,
                                     calibration*templates[1.0].image.array, rtol=1e-6)
        self.assertFloatsAlmostEqual(templates[calibration].variance.array,
                                     calibration**2*templates[1.0].variance.array, rtol=1e-6)

    def _remakePatches(self, photoCalib):
        """Rebuild the coadd patches from the main exposure with a different
        photometric calibration.
        """
        self.exposure.setPhotoCalib(photoCalib)
        self.patches.clear()
        self.dataIds.clear()
        for tract_id in range(4):
            self._makePatches(self.skymap.generateTract(tract_id))

    def testGaussianPsfGridFallback(self):
        """If the CoaddPsf cannot be evaluated at its average position, the
        Gaussian width is averaged over a grid of positions instead.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask()
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
        task = lsst.ip.diffim.GetTemplateTask()
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
        task = lsst.ip.diffim.GetTemplateTask()
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
        task = lsst.ip.diffim.GetTemplateTask()
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

    def _runCorrection(self, doCorrectVariancePlateScale=False, doScaleVariance=False,
                       handles=None):
        """Run the task on tract 0 with the named variance corrections.

        Everything defaults to off, so each test turns on exactly what it is
        exercising.
        """
        config = lsst.ip.diffim.GetTemplateTask.ConfigClass()
        config.doCorrectVariancePlateScale = doCorrectVariancePlateScale
        config.doScaleVariance = doScaleVariance
        task = lsst.ip.diffim.GetTemplateTask(config=config)
        result = task.run(coaddExposureHandles={0: handles or self.patches[0]},
                          bbox=lsst.geom.Box2I(self.varianceBox),
                          wcs=self.exposure.wcs,
                          dataIds={0: self.dataIds[0]},
                          physical_filter="a_test")
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
    def _makeCoaddInputs(records):
        """Make a CoaddInputs holding the given input records.

        Parameters
        ----------
        records : `list` [`tuple` [`lsst.afw.geom.SkyWcs` or `None`, \
                                   `lsst.geom.Box2I`]]
            The WCS and bbox to record for each input. A `None` WCS makes a
            record that `_plateScaleFactor` has to skip. May be empty, to
            simulate a coadd whose plate scale cannot be reconstructed.
        """
        ccdSchema = lsst.afw.table.ExposureTable.makeMinimalSchema()
        weightKey = ccdSchema.addField("weight", type=float, doc="Coadd weight")
        coaddInputs = lsst.afw.image.CoaddInputs(
            lsst.afw.table.ExposureTable.makeMinimalSchema(), ccdSchema)
        for wcs, bbox in records:
            record = coaddInputs.ccds.addNew()
            record.setWcs(wcs)
            record.setBBox(bbox)
            # Included because real coadds have it, though a single-record
            # correction does not use it.
            record.set(weightKey, 1.0)
        return coaddInputs

    def _patchHandles(self, tract, records):
        """Return handles for a tract's patches, with their CoaddInputs
        replaced by ``records``.
        """
        handles = []
        for ref in self.patches[tract]:
            coadd = ref.get()
            coadd.getInfo().setCoaddInputs(self._makeCoaddInputs(records))
            handles.append(pipeBase.InMemoryDatasetHandle(
                coadd, storageClass="ExposureF", copy=True, dataId=ref.dataId))
        return handles

    def testCorrectVariancePlateScale(self):
        """The plate scale correction is the total pixel area change from the
        images the coadds were built from to the science image.
        """
        _, off = self._runCorrection()

        # The fixture's records are the science image itself, so there is no
        # net change in pixel area: the same-instrument case.
        task, on = self._runCorrection(doCorrectVariancePlateScale=True)
        self.assertFloatsAlmostEqual(task.metadata["variancePlateScaleFactor"], 1.0, rtol=1e-6)

        # Coarser original pixels than science pixels, the DECam-template
        # case: the correction is the ratio of their areas.
        for scaleFactor in (1.315, 0.5):
            with self.subTest(scaleFactor=scaleFactor):
                handles = self._patchHandles(
                    0, [(self._makeScaledWcs(scaleFactor), self.exposure.getBBox())])
                task, on = self._runCorrection(doCorrectVariancePlateScale=True,
                                               handles=handles)
                factor = task.metadata["variancePlateScaleFactor"]
                self.assertFloatsAlmostEqual(factor, scaleFactor**2, rtol=1e-6)
                self.assertFloatsAlmostEqual(on.variance.array, off.variance.array*factor,
                                             rtol=1e-5, ignoreNaNs=True)

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
        task, _ = self._runCorrection(doCorrectVariancePlateScale=True, handles=handles)

        # The first record with a WCS: the one without is skipped, and the
        # last (which would give 1.0) is never reached.
        self.assertFloatsAlmostEqual(task.metadata["variancePlateScaleFactor"], 4.0, rtol=1e-6)

    def testCorrectVariancePlateScaleNeedsCoaddInputs(self):
        """Without usable coadd inputs the plate scale cannot be
        reconstructed, and the task must say so rather than silently applying
        only part of the correction.
        """
        handles = self._patchHandles(0, [])
        with self.assertRaisesRegex(RuntimeError, "doCorrectVariancePlateScale"):
            self._runCorrection(doCorrectVariancePlateScale=True, handles=handles)

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
            config = lsst.ip.diffim.GetTemplateTask.ConfigClass()
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
        taskOff, templateOff = _configureAndRunTask(False)
        self.assertFalse(hasattr(taskOff, "scaleVariance"))
        self.assertNotIn("scaleTemplateVarianceFactor", taskOff.metadata)

        # Both warps -- lanczos5 in ``_makePatches`` and lanczos3 in the
        # task -- correlate the noise. The variance plane tracks only the
        # per-pixel diagonal, which the second warp leaves too low, so
        # ``scaleVariance`` measures a factor well above 1 even though the
        # input variance planes are correct.
        #
        taskOn, templateOn = _configureAndRunTask(True)
        factor = taskOn.metadata["scaleTemplateVarianceFactor"]
        # TODO DM-55879: this value is pinned on purpose. The lanczos warping
        # kernels introduce small correlations that artificially suppress the
        # image pixel stddev and inflate the variance scaling factor. This
        # should be changed to 1.0 after DM-55879 is merged.
        self.assertFloatsAlmostEqual(factor, 1.1465, atol=0.01,
                                     msg="Measured template variance scaling changed; see the"
                                         " comment above if the correlation correction landed.")
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
        taskLow, templateLow = _configureAndRunTask(True, varianceScale=1/scaleFactor)
        self.assertFloatsAlmostEqual(taskLow.metadata["scaleTemplateVarianceFactor"],
                                     factor*scaleFactor, rtol=1e-5)
        self.assertImagesAlmostEqual(templateLow.variance, templateOn.variance, rtol=1e-5)

    def _runLegacyForFuture(self, raiseOnUndefinedMaskMap=True):
        """Build a template from tract 0 with a task configured for the
        future output type, without converting it.

        The task is kept as ``self.futureTask``.

        Parameters
        ----------
        raiseOnUndefinedMaskMap : `bool`, optional
            Value of the task config field of the same name.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            The legacy output struct.
        box : `lsst.geom.Box2I`
            The bounding box the template was requested on, before the task
            grew it by the template border.
        """
        config = lsst.ip.diffim.GetTemplateTask.ConfigClass()
        config.output_image_type = "future"
        config.raiseOnUndefinedMaskMap = raiseOnUndefinedMaskMap
        task = lsst.ip.diffim.GetTemplateTask(config=config)
        self.visit = 9876
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        result = task.run(coaddExposureHandles={0: self.patches[0]},
                          bbox=box,
                          wcs=self.exposure.wcs,
                          dataIds={0: self.dataIds[0]},
                          physical_filter="a_test",
                          visit=self.visit)
        # run() always produces legacy types, whatever output_image_type is.
        self.assertIsInstance(result.template, lsst.afw.image.ExposureF)
        self.futureTask = task
        return result, box

    def _runFuture(self):
        """Build a template and convert it, the way ``runQuantum`` does.

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
        result, box = self._runLegacyForFuture()
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
                self.assertFalse(template.psf_shape_flag)
                self.assertGreater(template.psf_shape_xx, 0)
                self.assertGreater(template.bounds.area, 0)

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
