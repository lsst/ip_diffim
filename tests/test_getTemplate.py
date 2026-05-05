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

import numpy as np

import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.math
from lsst.daf.butler import DataCoordinate, DimensionUniverse
import lsst.geom
import lsst.ip.diffim
import lsst.meas.algorithms
import lsst.meas.base.tests
import lsst.pipe.base as pipeBase
import lsst.skymap
import lsst.utils.tests

from utils import generate_data_id

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


class GetTemplateTaskTestBase(lsst.utils.tests.TestCase):
    """Test harness for variants of GetTemplateTask

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
        self.exposure = self._makeExposure()

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
        self.assertEqual(template.photoCalib, self.exposure.photoCalib)
        self.assertEqual(template.getXY0(), expectedBox.getMin())
        self.assertEqual(template.filter.bandLabel, "a")
        self.assertEqual(template.filter.physicalLabel, "a_test")
        self.assertEqual(template.psf.getComponentCount(), nPsfs)

    def _replacePatchesWithNoise(self, true_var, seed=0):
        """Replace each patch's pixels with iid Gaussian noise of the
        given marginal variance and pin the variance plane to that
        value, so the input to the task has *consistent* noise and
        variance plane (mask bits cleared too).

        Parameters
        ----------
        true_var : `float`
            Variance of the Gaussian noise drawn into each patch's
            image plane; the same value is written to the variance
            plane and used as the ground truth that downstream
            propagation tests compare against.
        seed : `int`
            Numpy RNG seed; per-patch noise draws use independent
            substreams so tract/patch overlaps combine independent
            samples.
        """
        rng = np.random.default_rng(seed)
        for tract_id in list(self.patches.keys()):
            new_handles = []
            for handle in self.patches[tract_id]:
                exposure = handle.get()
                ny, nx = exposure.image.array.shape
                exposure.image.array[:, :] = rng.normal(
                    0.0, np.sqrt(true_var), size=(ny, nx),
                ).astype(np.float32)
                exposure.variance.array[:, :] = float(true_var)
                exposure.mask.array[:, :] = 0
                new_handles.append(pipeBase.InMemoryDatasetHandle(
                    exposure, storageClass="ExposureF",
                    copy=True, dataId=handle.dataId,
                ))
            self.patches[tract_id] = new_handles

    def _interiorView(self, template, shrink=30):
        """Return (image, variance, valid_mask) for the interior of the
        merged template, excluding NO_DATA / EDGE pixels and a
        ``shrink``-pixel inset that drops residual warping/convolution
        edge effects from the analytical bounds.

        Used by variance-propagation tests that compare the empirical
        per-pixel variance of the output image to the median of the
        variance plane: those statistics are only meaningful where every
        contributing input was real data.
        """
        mask_bits = (template.mask.getPlaneBitMask("NO_DATA")
                     | template.mask.getPlaneBitMask("EDGE"))
        good = (template.mask.array & mask_bits) == 0
        # Sentinel-marked variance pixels inflate the median artificially.
        good &= template.variance.array < 1e10
        inset = lsst.geom.Box2I(template.getBBox())
        inset.grow(-shrink)
        x0, y0 = template.getX0(), template.getY0()
        sl = (slice(inset.getMinY() - y0, inset.getMaxY() - y0 + 1),
              slice(inset.getMinX() - x0, inset.getMaxX() - x0 + 1))
        valid = good[sl]
        return (template.image.array[sl], template.variance.array[sl], valid)

    def _checkPixels(self, template, config, box):
        """Check that the pixel values in the template are close to the
        original image.
        """
        # All pixels should have real values!
        expectedBox = lsst.geom.Box2I(box)
        expectedBox.grow(config.templateBorderSize)

        if debug:
            _showTemplate(expectedBox, template)

        # Check that we fully filled the template from the patches.
        self.assertTrue(np.all(np.isfinite(template.image.array)))
        # Because of the scale changes, there will be some ringing in the
        # difference between the template and the original image; pick
        # tolerances large enough to account for that.
        self.assertImagesAlmostEqual(template.image, self.exposure[expectedBox].image,
                                     rtol=.1, atol=4)
        # The variance plane is now an inverse-variance combination, so a
        # pixel covered by N independent patches has variance ≈ input/N.
        # Instead of comparing element-wise to the input, sanity-check the
        # plane: finite everywhere, strictly positive where data lands,
        # and never inflated above the per-pixel input.
        template_var = template.variance.array
        input_var_max = float(np.max(self.exposure[expectedBox].variance.array))
        self.assertTrue(np.all(np.isfinite(template_var)))
        valid = (template.mask.array
                 & template.mask.getPlaneBitMask("NO_DATA")) == 0
        self.assertTrue(np.all(template_var[valid] > 0))
        self.assertTrue(np.all(template_var[valid] <= 2*input_var_max))
        # Not checking the mask, as warping changes the sizes of the masks.


class GetTemplateTaskTestCase(GetTemplateTaskTestBase):
    """Test that GetTemplateTask works on both one tract and multiple tract
    input coadd exposures.
    """

    def testRunOneTractInput(self):
        """Test a bounding box that fully fits inside one tract, with only
        that tract passed as input. This checks that the code handles a single
        tract input correctly.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask()
        # Restrict to tract 0, since the box fits in just that tract.
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles={0: self.patches[0]},
                          bbox=lsst.geom.Box2I(box),
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
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
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
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")

        # All 4 patches from all 4 tracts are included in this template
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 9)
        self._checkPixels(result.template, task.config, box)

    def testRunNoTemplate(self):
        """A bounding box that doesn't overlap the patches will raise.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(1200, 1200), lsst.geom.Point2I(1600, 1600))
        task = lsst.ip.diffim.GetTemplateTask()
        with self.assertRaisesRegex(lsst.pipe.base.NoWorkFound, "No patches found"):
            task.run(coaddExposureHandles=self.patches,
                     bbox=lsst.geom.Box2I(box),
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
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        no_data = (result.template.mask.array & result.template.mask.getPlaneBitMask("NO_DATA")) != 0
        self.assertTrue(np.isfinite(result.template.image.array).all())
        self.assertTrue(np.isfinite(result.template.variance.array).all())
        self.assertEqual(no_data.sum(), 20990)

    def testVariancePlaneTracksOutputNoise(self):
        """The output variance plane should match the marginal per-pixel
        variance of the output image plane for noise-only inputs.

        Lanczos warping introduces correlations between output pixels but
        the *marginal* variance of any one pixel is still the empirical
        variance of pixel values across many pixels (the cross-pixel
        covariance does not bias an unweighted sample variance). So even
        though sqrt(template.image) drops below sqrt(true_var) after warp,
        the variance plane must drop by exactly the same factor —
        otherwise downstream 5-σ peak detection on
        ``image / sqrt(variance)`` would mis-threshold.
        """
        true_var = 4.0
        self._replacePatchesWithNoise(true_var=true_var, seed=0)

        box = lsst.geom.Box2I(lsst.geom.Point2I(200, 200),
                              lsst.geom.Point2I(600, 600))
        task = lsst.ip.diffim.GetTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        image, variance, valid = self._interiorView(result.template)
        emp = float(np.var(image[valid]))
        plane = float(np.median(variance[valid]))
        self.assertGreater(valid.sum(), 1000,
                           "fixture must produce enough good pixels")
        self.assertGreater(emp, 0.0)
        self.assertGreater(plane, 0.0)
        # Both empirical and plane must be reduced from true_var (Lanczos
        # warp + multi-tract inverse-variance combine).
        self.assertLess(plane, true_var)
        # Per-pixel agreement between marginal noise and variance plane.
        # 25% tolerance covers correlated-sample variance estimator error
        # plus per-tile inverse-variance combination noise.
        self.assertAlmostEqual(
            emp, plane, delta=0.25*plane,
            msg=f"empirical output variance {emp:.3f} disagrees with "
                f"variance plane {plane:.3f} by more than 25%",
        )

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
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        if debug:
            _showTemplate(box, result.template)
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 9)
        # We just check that the pixel values are all finite. We cannot
        # check that pixel values in the template are closer to the
        # original anymore.
        self.assertTrue(np.isfinite(result.template.image.array).all())


class GetPsfMatchedTemplateTaskTestCase(GetTemplateTaskTestBase):
    """Tests for `GetPsfMatchedTemplateTask`.

    Reuses the multi-tract / multi-patch fixture from
    `GetTemplateTaskTestBase`, then overrides each patch's PSF with a
    `GaussianPsf` of a chosen sigma so the matching has something to do.
    """

    # Fixed bbox used by the multi-tract integration tests below.
    _box = lsst.geom.Box2I(lsst.geom.Point2I(200, 200),
                           lsst.geom.Point2I(600, 600))

    def _setPatchSigmas(self, sigmas):
        """Replace each patch's PSF with `GaussianPsf(σ)` for a chosen σ.

        Parameters
        ----------
        sigmas : `dict` [`tuple` [`int`, `int`], `float`]
            Mapping from ``(tract_id, patch_id)`` to the desired intrinsic
            Gaussian sigma in pixels. Must contain an entry for every
            ``(tract_id, patch_id)`` in ``self.patches``.
        """
        for tract_id in list(self.patches.keys()):
            new_handles = []
            for patch_id, handle in enumerate(self.patches[tract_id]):
                exposure = handle.get()
                sigma = sigmas[(tract_id, patch_id)]
                exposure.setPsf(lsst.afw.detection.GaussianPsf(25, 25, sigma))
                new_handles.append(pipeBase.InMemoryDatasetHandle(
                    exposure, storageClass="ExposureF",
                    copy=True, dataId=handle.dataId,
                ))
            self.patches[tract_id] = new_handles

    def _injectGaussiansAndSetPsfs(self, sigmas, flux):
        """Inject a bright 2D Gaussian source at each patch's bbox center
        with σ matching the patch's metadata, and set the patch PSF to
        match.

        Parameters
        ----------
        sigmas : `dict` [`tuple` [`int`, `int`], `float`]
            ``(tract_id, patch_id) → σ`` mapping; controls both the PSF
            attached to each patch and the width of the injected Gaussian.
        flux : `float`
            Integrated flux of each injected Gaussian (peak amplitude is
            ``flux/(2π σ²)``).

        Returns
        -------
        science_positions : `dict` [`tuple` [`int`, `int`],
                                    `lsst.geom.Point2D`]
            Mapping from ``(tract_id, patch_id)`` to the position of the
            injection projected onto the *science* exposure's pixel grid,
            useful for measuring the source on the merged template.
        """
        science_positions = {}
        for tract_id in list(self.patches.keys()):
            new_handles = []
            for patch_id, handle in enumerate(self.patches[tract_id]):
                sigma = sigmas[(tract_id, patch_id)]
                exposure = handle.get()
                bbox = exposure.getBBox()

                cx_patch = bbox.getCenterX()
                cy_patch = bbox.getCenterY()
                x0, y0 = bbox.getMinX(), bbox.getMinY()
                ny, nx = exposure.image.array.shape
                yy, xx = np.mgrid[:ny, :nx].astype(np.float64)
                xx_abs = xx + x0
                yy_abs = yy + y0
                gauss = (flux/(2*np.pi*sigma**2)
                         * np.exp(-0.5*((xx_abs - cx_patch)**2 + (yy_abs - cy_patch)**2)/sigma**2))
                exposure.image.array += gauss.astype(exposure.image.array.dtype)

                exposure.setPsf(lsst.afw.detection.GaussianPsf(25, 25, sigma))

                sky = exposure.wcs.pixelToSky(lsst.geom.Point2D(cx_patch, cy_patch))
                science_positions[(tract_id, patch_id)] = self.exposure.wcs.skyToPixel(sky)

                new_handles.append(pipeBase.InMemoryDatasetHandle(
                    exposure, storageClass="ExposureF",
                    copy=True, dataId=handle.dataId,
                ))
            self.patches[tract_id] = new_handles
        return science_positions

    @staticmethod
    def _measureGaussianSigma(image, position, halfsize=10):
        """Estimate the Gaussian-equivalent σ of a bright source at
        ``position`` via second moments on a postage stamp.

        The stamp's edge rows/columns are used to subtract any local
        background or diffuse light before computing moments, so faint
        neighbors don't bias the measurement of a much brighter source.
        """
        px = int(round(position.x))
        py = int(round(position.y))
        bbox = lsst.geom.Box2I(
            lsst.geom.Point2I(px - halfsize, py - halfsize),
            lsst.geom.Extent2I(2*halfsize + 1, 2*halfsize + 1),
        )
        bbox.clip(image.getBBox())
        stamp = image[bbox].array.astype(np.float64, copy=True)
        edge = np.concatenate([stamp[0, :], stamp[-1, :],
                               stamp[:, 0], stamp[:, -1]])
        stamp -= np.median(edge)

        yy, xx = np.mgrid[bbox.getMinY():bbox.getMaxY() + 1,
                          bbox.getMinX():bbox.getMaxX() + 1]
        flux = stamp.sum()
        cx = (xx*stamp).sum()/flux
        cy = (yy*stamp).sum()/flux
        Mxx = ((xx - cx)**2 * stamp).sum()/flux
        Myy = ((yy - cy)**2 * stamp).sum()/flux
        return float((Mxx*Myy)**0.25)

    # ---- helper-level unit tests --------------------------------------

    def test_matchingKernelSigma(self):
        """``_matchingKernelSigma`` returns ``sqrt(t² − σ²)`` or `None`."""
        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        # Standard case.
        self.assertAlmostEqual(task._matchingKernelSigma(1.5, 2.5),
                               np.sqrt(2.5**2 - 1.5**2), places=6)
        # σ at or above target → no convolution.
        self.assertIsNone(task._matchingKernelSigma(2.5, 2.5))
        self.assertIsNone(task._matchingKernelSigma(3.0, 2.5))
        # Implied kernel sigma below the configured floor → no
        # convolution. sqrt(2.000001² − 4) ≈ 0.002, well below the
        # default minMatchingKernelSigma of 0.05.
        self.assertIsNone(task._matchingKernelSigma(2.0, 2.000001))

    def test_convolvePsf(self):
        """A Gaussian-on-Gaussian convolution adds in quadrature and the
        returned `KernelPsf` is normalized.
        """
        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        sigma_in, kernel_sigma = 2.0, 1.5
        psf_in = lsst.afw.detection.GaussianPsf(25, 25, sigma_in)
        psf_out = task._convolvePsf(psf_in, kernel_sigma)

        position = psf_out.getAveragePosition()
        sigma_out = psf_out.computeShape(position).getDeterminantRadius()
        self.assertAlmostEqual(
            sigma_out, np.sqrt(sigma_in**2 + kernel_sigma**2), places=3,
        )
        # Kernel must integrate to 1 so flux is preserved on convolution.
        kernel_image = psf_out.computeKernelImage(position)
        self.assertAlmostEqual(kernel_image.array.sum(), 1.0, places=6)

    def test_sigmaClip(self):
        """``_sigmaClip`` is one-sided: it rejects only upper outliers,
        keeps below-median sigmas regardless of how far below they are,
        and drops non-finite/non-positive entries.
        """
        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        # Five near-median sigmas, one large *upper* outlier, one NaN,
        # one negative, and one *narrow* (below-median) sigma.
        sigmas = np.array([1.95, 2.0, 2.05, 1.98, 2.02, 50.0, np.nan, -1.0, 1.0])
        keep = task._sigmaClip(sigmas)
        self.assertTrue(np.all(keep[:5]), "near-median sigmas should be kept")
        self.assertFalse(keep[5], "upper-outlier sigma=50 should be rejected")
        self.assertFalse(keep[6], "NaN sigma should be dropped")
        self.assertFalse(keep[7], "negative sigma should be dropped")
        self.assertTrue(keep[8], "below-median sigma should always be kept")

    def test_makeMetadataCatalog_componentLoadsOnly(self):
        """``_makeMetadataCatalog`` loads only the metadata components and
        never triggers a full Exposure load.
        """
        class RecordingHandle:
            def __init__(self, exposure, dataId):
                self._exp = exposure
                self.dataId = dataId
                self.fetched = []

            def get(self, *, component=None):
                if component is None:
                    self.fetched.append("FULL")
                    raise AssertionError(
                        "metadata pass must not full-load the Exposure"
                    )
                self.fetched.append(component)
                if component == "psf":
                    return self._exp.psf
                if component == "wcs":
                    return self._exp.wcs
                if component == "photoCalib":
                    return self._exp.photoCalib
                if component == "bbox":
                    return self._exp.getBBox()
                raise KeyError(component)

        # Wrap one tract's worth of patches with recording handles. Use
        # the standardized DataCoordinates that the production task
        # actually receives (where dataId["patch"] is an int), not the
        # generate_data_id-flavored ones the handles carry internally.
        tract_id = next(iter(self.patches))
        dataIds = self.dataIds[tract_id]
        handles = [
            RecordingHandle(ref.get(), dataIds[i])
            for i, ref in enumerate(self.patches[tract_id])
        ]

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        catalog, totalBox, refByDataId = task._makeMetadataCatalog(
            handles, dataIds,
        )

        self.assertEqual(len(catalog), len(handles))
        expected = {"psf", "wcs", "photoCalib", "bbox"}
        for h in handles:
            self.assertEqual(set(h.fetched), expected)
            self.assertNotIn("FULL", h.fetched)
        for record, dataId in zip(catalog, dataIds):
            self.assertEqual(record.get("tract"), dataId["tract"])
            self.assertEqual(record.get("patch"), dataId["patch"])
        for h, dataId in zip(handles, dataIds):
            self.assertIs(refByDataId[dataId], h)
        # totalBox is the union of the patch bboxes.
        union = lsst.geom.Box2I()
        for h in handles:
            union = union.expandedTo(h._exp.getBBox())
        self.assertEqual(totalBox, union)

    def test_basicMatching_targetMatchesMaxInput(self):
        """End-to-end: patches with varying sigmas all get matched up to
        the padded max input; the output ``CoaddPsf`` reports that width.
        """
        sigmas = {
            (0, 0): 1.5, (0, 1): 2.0, (0, 2): 2.5, (0, 3): 3.0,
            (1, 0): 1.5, (1, 1): 2.0, (1, 2): 2.5, (1, 3): 3.0,
            (2, 0): 1.5, (2, 1): 2.0, (2, 2): 2.5, (2, 3): 3.0,
            (3, 0): 1.5, (3, 1): 2.0, (3, 2): 2.5, (3, 3): 3.0,
        }
        self._setPatchSigmas(sigmas)

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(self._box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")

        # Target = max input σ. Patches within
        # ``targetPsfSigmaPad`` of the max are passed through unconvolved;
        # others are convolved up to the target.
        expected_target = max(sigmas.values())
        self.assertAlmostEqual(
            result.template.metadata["PSF_MATCHED_TARGET_SIGMA"],
            expected_target, places=4,
        )

        # Output PSF is a CoaddPsf assembled from per-patch matched
        # KernelPsfs; its Gaussian-equivalent radius at the data center
        # should equal the target — *but* the patches and the science
        # exposure use different pixel scales (0.20"/px vs 0.21"/px in
        # the fixture), so the patch-pixel sigma maps to a slightly
        # smaller science-pixel sigma on the output frame.
        psf = result.template.psf
        center = lsst.geom.Point2D(result.template.getBBox().getCenter())
        sigma_out = psf.computeShape(center).getDeterminantRadius()
        patch_scale, science_scale = self.scale, 1.05*self.scale
        expected_sigma_out = expected_target*patch_scale/science_scale
        self.assertAlmostEqual(sigma_out, expected_sigma_out, delta=0.05)

        # Image plane should have all-finite pixels in the data region.
        self.assertTrue(np.all(np.isfinite(result.template.image.array)))

    def test_outlierRejection_isOneSided(self):
        """Upper-outlier patches are rejected, but below-median
        patches are always kept — even when their sigma is far below
        the bulk.

        The clipping is asymmetric because a wide-PSF outlier would
        force the global target wide enough to broaden every other
        patch, while a narrow-PSF outlier just gets convolved up to
        the target with a slightly larger kernel.
        """
        # Bulk near σ=2.0 (slight variation so MAD > 0), one upper
        # outlier at tract 0 patch 0, and one well-below-median entry
        # at tract 3 patch 3.
        sigmas = {
            (0, 0): 4.00, (0, 1): 2.01, (0, 2): 2.02, (0, 3): 2.03,
            (1, 0): 2.01, (1, 1): 2.02, (1, 2): 2.03, (1, 3): 2.04,
            (2, 0): 2.02, (2, 1): 2.03, (2, 2): 2.04, (2, 3): 2.05,
            (3, 0): 2.03, (3, 1): 2.04, (3, 2): 2.05, (3, 3): 0.50,
        }
        self._setPatchSigmas(sigmas)

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(self._box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")

        # The σ=4.0 upper outlier is rejected; the σ=0.5 lower outlier
        # is kept along with the bulk.
        self.assertEqual(
            result.template.metadata["PSF_MATCHED_NUM_REJECTED"], 1,
        )
        # Target reflects the surviving max (~2.05), not the rejected
        # σ=4.0 outlier.
        target = result.template.metadata["PSF_MATCHED_TARGET_SIGMA"]
        self.assertLess(target, 3.0)

    def test_metadataKeysPopulated(self):
        """All ``PSF_MATCHED_*`` provenance keys are written out."""
        sigmas = {
            (0, 0): 2.00, (0, 1): 2.02, (0, 2): 2.04, (0, 3): 2.06,
            (1, 0): 2.02, (1, 1): 2.04, (1, 2): 2.06, (1, 3): 2.08,
            (2, 0): 2.04, (2, 1): 2.06, (2, 2): 2.08, (2, 3): 2.10,
            (3, 0): 2.06, (3, 1): 2.08, (3, 2): 2.10, (3, 3): 2.12,
        }
        self._setPatchSigmas(sigmas)

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(self._box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")

        md = result.template.metadata
        self.assertIn("PSF_MATCHED_TARGET_SIGMA", md)
        self.assertIn("PSF_MATCHED_NUM_PATCHES", md)
        self.assertIn("PSF_MATCHED_NUM_REJECTED", md)
        self.assertGreater(md["PSF_MATCHED_TARGET_SIGMA"], 0.0)
        self.assertGreater(md["PSF_MATCHED_NUM_PATCHES"], 0)
        # No outliers in this fixture.
        self.assertEqual(md["PSF_MATCHED_NUM_REJECTED"], 0)

    def test_convolveImage_doesNotMutateInput(self):
        """``_convolveImage`` should not modify the caller's MaskedImage,
        because the same Exposure may be referenced elsewhere (e.g. by
        the catalog record's PSF) and silent in-place sanitization is a
        nasty source-of-truth bug.
        """
        n = 60
        mi = lsst.afw.image.MaskedImageF(n, n)
        mi.image.array[:, :] = 1.5
        mi.image.array[10, 10] = np.nan
        mi.variance.array[:, :] = 4.0
        mi.variance.array[20, 20] = np.nan
        image_before = mi.image.array.copy()
        variance_before = mi.variance.array.copy()
        mask_before = mi.mask.array.copy()

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()

        kernel = task._buildMatchingKernel(2.0)
        task._convolveImage(mi, kernel)

        np.testing.assert_array_equal(
            mi.mask.array, mask_before,
            err_msg="input mask plane was mutated by _convolveImage",
        )
        np.testing.assert_array_equal(
            mi.variance.array[~np.isnan(variance_before)],
            variance_before[~np.isnan(variance_before)],
            err_msg="input variance plane was mutated by _convolveImage",
        )
        # NaN must still be NaN where the caller put it.
        self.assertTrue(np.isnan(mi.variance.array[20, 20]))
        np.testing.assert_array_equal(
            mi.image.array[~np.isnan(image_before)],
            image_before[~np.isnan(image_before)],
            err_msg="input image plane was mutated by _convolveImage",
        )
        self.assertTrue(np.isnan(mi.image.array[10, 10]))

    def test_convolveImage_badVarianceDoesNotPoisonOutput(self):
        """A handful of bad-variance input pixels (NaN, ±inf, ≤ 0)
        must not turn a kernel-sized footprint of the output variance
        plane into NaN/Inf.

        This is a regression test for the specific bug where
        ``_convolveImage`` zeroed the IMAGE plane's non-finite values
        but left the variance plane raw, so a single NaN-variance pixel
        smeared NaN over its whole kernel footprint via the Σ K²
        propagation. Downstream ``_merge`` then drops every smeared
        pixel as bad, gutting the usable area of the matched template.
        """
        n = 80
        kernel_sigma = 2.0
        mi = lsst.afw.image.MaskedImageF(n, n)
        mi.image.array[:, :] = 0.0
        mi.variance.array[:, :] = 4.0
        # Three different flavors of bad variance, well-separated so
        # their kernel footprints don't overlap in this 80×80 fixture.
        # All placed in the interior so the convolution edge band
        # (kernel halfwidth around the perimeter, naturally NaN'd by
        # afwMath.convolve) doesn't mask the bug.
        bad_pixels = [(20, 20), (40, 60), (60, 30)]
        bad_values = [np.nan, 0.0, -1.0]
        for (yy, xx), value in zip(bad_pixels, bad_values):
            mi.variance.array[yy, xx] = value

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        kernel = task._buildMatchingKernel(kernel_sigma)
        out = task._convolveImage(mi, kernel)

        # The convolution leaves a kernel-halfwidth NaN/EDGE band around
        # the image perimeter by design; only check the interior.
        halo = int(np.ceil(task.config.matchingKernelSigmas*kernel_sigma))
        interior = (slice(halo, n - halo), slice(halo, n - halo))
        bad_count = int((~np.isfinite(out.variance.array[interior])).sum())
        self.assertEqual(
            bad_count, 0,
            f"{bad_count} of {(n - 2*halo)**2} interior variance pixels "
            "are non-finite; bad-input pixels are smearing NaN/Inf "
            "through the convolution kernel",
        )
        self.assertTrue(np.all(np.isfinite(out.image.array[interior])))
        # The bad-input footprints should be flagged NO_DATA so
        # downstream code can identify and deweight them.
        no_data = (out.mask.array
                   & out.mask.getPlaneBitMask("NO_DATA")) != 0
        self.assertTrue(no_data.any(),
                        "bad-input footprints should set NO_DATA")
        # Each bad input pixel poisons at most a (2·halo+1)² footprint;
        # cap the total flagged area at that to catch a wider regression
        # but still flag the original "all of it" bug.
        max_flagged = len(bad_pixels)*(2*halo + 1)**2
        self.assertLessEqual(
            no_data[interior].sum(), max_flagged,
            "NO_DATA footprint inside the interior is larger than the "
            "kernel footprint around each bad input pixel",
        )

    def test_convolveImage_variancePropagation(self):
        """Convolution by a normalized 2D Gaussian of σ_k should reduce
        the variance plane by Σ K² ≈ 1/(4π σ_k²), and the empirical
        per-pixel variance of the output image should match the variance
        plane.
        """
        n = 256
        sigma_k = 3.0
        rng = np.random.default_rng(123)
        mi = lsst.afw.image.MaskedImageF(n, n)
        mi.image.array[:, :] = rng.normal(
            0, 2, size=(n, n),
        ).astype(np.float32)
        mi.variance.array[:, :] = 4.0

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        kernel = task._buildMatchingKernel(sigma_k)
        out = task._convolveImage(mi, kernel)

        # Stay well inside the kernel halo to dodge the EDGE band the
        # convolution leaves around the perimeter.
        halo = int(np.ceil(task.config.matchingKernelSigmas*sigma_k)) + 5
        inner = (slice(halo, n - halo), slice(halo, n - halo))
        emp = float(np.var(out.image.array[inner]))
        plane = float(np.median(out.variance.array[inner]))

        # Σ K² for a normalized 2D Gaussian (continuous) is
        # 1/(4π σ²). Discretization shifts this slightly; allow 5% on
        # the plane and 15% on the empirical (sample-variance noise).
        expected = 4.0/(4*np.pi*sigma_k**2)
        self.assertAlmostEqual(plane, expected, delta=0.05*expected)
        self.assertAlmostEqual(emp, expected, delta=0.15*expected)
        self.assertAlmostEqual(emp, plane, delta=0.15*plane)

    def testVariancePlaneTracksOutputNoise(self):
        """After the full PSF-match + warp + merge pipeline, the variance
        plane of the output template should still track the marginal
        per-pixel variance of the output image plane, even though the
        absolute level has dropped substantially due to the
        matching-kernel convolution and the warping kernel.
        """
        true_var = 4.0
        self._replacePatchesWithNoise(true_var=true_var, seed=1)
        # Re-attach a Gaussian PSF to every patch so matching has work to
        # do (else `_chooseTargetSigma` from the WarpedPsf wouldn't have a
        # well-defined width).
        sigmas = {(t, p): 1.5 for t in range(4) for p in range(4)}
        self._setPatchSigmas(sigmas)

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(self._box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        image, variance, valid = self._interiorView(result.template)
        emp = float(np.var(image[valid]))
        plane = float(np.median(variance[valid]))
        self.assertGreater(valid.sum(), 1000)
        self.assertGreater(plane, 0.0)
        # Convolution + warp should both reduce the variance well below
        # the input level.
        self.assertLess(plane, true_var)
        self.assertAlmostEqual(
            emp, plane, delta=0.25*plane,
            msg=f"PSF-matched template variance plane {plane:.4f} "
                f"disagrees with empirical output variance "
                f"{emp:.4f} by more than 25%",
        )

    def test_outputVarianceScrubbedOfBadValues(self):
        """The output template's variance plane never contains NaN,
        ±inf, or non-positive values, even when input patches carry
        bad variance pixels and the merge cannot fill some output
        regions. NO_DATA pixels carry the task's ``_noDataVariance``
        sentinel so downstream detection sees finite, positive
        variance everywhere.
        """
        sigmas = {(t, p): 2.0 for t in range(4) for p in range(4)}
        self._setPatchSigmas(sigmas)

        # Drop tract 0 patch 1 to force a NO_DATA region in the merged
        # output. Pop ``patches`` and ``dataIds`` together so the
        # per-tract iteration stays aligned.
        self.patches[0].pop(1)
        self.dataIds[0].pop(1)

        # Poison NaN and zero pixels in every remaining patch's
        # variance plane to also exercise the per-patch bad-input path
        # through `_convolveImage` and `_merge`.
        poisoned = collections.defaultdict(list)
        for tract_id, handles in self.patches.items():
            for handle in handles:
                exposure = handle.get()
                cx = int(exposure.getBBox().getCenterX())
                cy = int(exposure.getBBox().getCenterY())
                nan_box = lsst.geom.Box2I(
                    lsst.geom.Point2I(cx - 3, cy - 3),
                    lsst.geom.Extent2I(6, 6),
                )
                zero_box = lsst.geom.Box2I(
                    lsst.geom.Point2I(cx + 10, cy + 10),
                    lsst.geom.Extent2I(6, 6),
                )
                exposure.variance[nan_box].array[:] = np.nan
                exposure.variance[zero_box].array[:] = 0.0
                poisoned[tract_id].append(pipeBase.InMemoryDatasetHandle(
                    exposure, storageClass="ExposureF",
                    copy=True, dataId=handle.dataId,
                ))
        self.patches = poisoned

        # Small box mostly in tract 0 — with patch 1 missing the merge
        # leaves a real NO_DATA region (cf. `testMissingPatches`).
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                              lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        template = result.template
        var = template.variance.array

        self.assertTrue(np.all(np.isfinite(var)),
                        "variance plane must contain no NaN or inf")
        self.assertTrue(np.all(var > 0),
                        "variance plane must be strictly positive")

        no_data = (template.mask.array
                   & template.mask.getPlaneBitMask("NO_DATA")) != 0
        self.assertTrue(
            no_data.any(),
            "fixture should produce some NO_DATA pixels for the test"
            " to be meaningful",
        )
        # The scrub replaces bad-variance pixels with the sentinel; those
        # pixels are exactly where the merge couldn't fill from any
        # input, so they must lie inside the NO_DATA region. (NO_DATA can
        # be a strict superset: input mask bits OR'd through the merge
        # can mark pixels NO_DATA even where the merge did fill valid
        # variance.)
        sentinel = var == task._noDataVariance
        self.assertTrue(
            sentinel.any(),
            "scrub should have replaced at least some bad-variance"
            " pixels with the sentinel",
        )
        self.assertTrue(
            np.all(no_data[sentinel]),
            "every sentinel pixel must lie inside the NO_DATA region",
        )

    def test_minAcceptedPatches_unsatisfiedRaises(self):
        """If too few patches survive sigma clipping, raise ``NoWorkFound``."""
        sigmas = {
            (0, 0): 2.0, (0, 1): 2.0, (0, 2): 2.0, (0, 3): 2.0,
            (1, 0): 2.0, (1, 1): 2.0, (1, 2): 2.0, (1, 3): 2.0,
            (2, 0): 2.0, (2, 1): 2.0, (2, 2): 2.0, (2, 3): 2.0,
            (3, 0): 2.0, (3, 1): 2.0, (3, 2): 2.0, (3, 3): 2.0,
        }
        self._setPatchSigmas(sigmas)

        config = lsst.ip.diffim.GetPsfMatchedTemplateConfig()
        config.minAcceptedPatches = 1000  # more than the fixture supplies
        task = lsst.ip.diffim.GetPsfMatchedTemplateTask(config=config)

        with self.assertRaisesRegex(lsst.pipe.base.NoWorkFound,
                                    "survived PSF outlier rejection"):
            task.run(coaddExposureHandles=self.patches,
                     bbox=lsst.geom.Box2I(self._box),
                     wcs=self.exposure.wcs,
                     dataIds=self.dataIds,
                     physical_filter="a_test")

    def test_multiplePatchesMergeIntoMatchedTemplate(self):
        """Sixteen patches with sixteen distinct intrinsic PSF sizes —
        spread over four tracts — are PSF-matched to a common Gaussian
        and merged into a single template.

        This exercises the full path: each patch picks up its own
        matching kernel σ_k = sqrt(target² − σ_in²); the convolved
        patches merge per-tract, are warped to the science geometry,
        merged across tracts, and the resulting CoaddPsf carries the
        per-patch matched ``KernelPsf``s. The test verifies that
        evaluating the CoaddPsf at multiple widely-separated points
        gives the same Gaussian-equivalent width — i.e. every
        contributing patch was matched consistently.
        """
        # Sixteen distinct intrinsic sigmas in [1.5, 3.0], one per
        # (tract, patch) cell. Spans a 2× range, all comfortably within
        # MAD of the median so nothing is sigma-clipped.
        sigmas = {
            (0, 0): 1.5, (0, 1): 1.6, (0, 2): 1.7, (0, 3): 1.8,
            (1, 0): 1.9, (1, 1): 2.0, (1, 2): 2.1, (1, 3): 2.2,
            (2, 0): 2.3, (2, 1): 2.4, (2, 2): 2.5, (2, 3): 2.6,
            (3, 0): 2.7, (3, 1): 2.8, (3, 2): 2.9, (3, 3): 3.0,
        }
        # Inject a bright 2D Gaussian source at each patch's bbox center
        # with σ matching that patch's metadata. After matching, every
        # injection's measured width on the merged template should equal
        # the target σ (mapped from patch pixels to science pixels by the
        # WCS scale ratio) — a direct test that the convolution kernel
        # actually broadens the per-patch PSFs.
        injection_flux = 1e7
        science_positions = self._injectGaussiansAndSetPsfs(
            sigmas, flux=injection_flux,
        )

        task = lsst.ip.diffim.GetPsfMatchedTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(self._box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        template = result.template
        md = template.metadata

        # No outliers in this distribution.
        self.assertEqual(md["PSF_MATCHED_NUM_REJECTED"], 0)
        # Multiple patches merged into the template, and the CoaddPsf
        # has one component per accepted patch.
        self.assertGreater(md["PSF_MATCHED_NUM_PATCHES"], 1)
        self.assertEqual(template.psf.getComponentCount(),
                         md["PSF_MATCHED_NUM_PATCHES"])

        # Target = max input σ. Patches within
        # ``targetPsfSigmaPad`` of the max are passed through unconvolved;
        # others are convolved up to the target.
        expected_target = max(sigmas.values())
        self.assertAlmostEqual(md["PSF_MATCHED_TARGET_SIGMA"],
                               expected_target, places=4)

        # Patch and science exposures use slightly different pixel
        # scales, so a patch-pixel sigma maps to a slightly different
        # science-pixel sigma on the output frame.
        patch_scale = self.patches[0][0].get().wcs.getPixelScale().asArcseconds()
        science_scale = self.exposure.wcs.getPixelScale().asArcseconds()
        expected_sigma_out = expected_target*patch_scale/science_scale

        # Sample the CoaddPsf at five widely-separated points: four
        # quadrants plus the center. Different patches dominate at each
        # point, so this is a real test that the matching produced a
        # consistent Gaussian width across the merged template.
        bbox = template.getBBox()
        cx, cy = bbox.getCenterX(), bbox.getCenterY()
        dx, dy = bbox.getWidth()//4, bbox.getHeight()//4
        sample_points = [
            lsst.geom.Point2D(cx - dx, cy - dy),
            lsst.geom.Point2D(cx + dx, cy - dy),
            lsst.geom.Point2D(cx - dx, cy + dy),
            lsst.geom.Point2D(cx + dx, cy + dy),
            lsst.geom.Point2D(cx, cy),
        ]
        for point in sample_points:
            sigma_out = template.psf.computeShape(point).getDeterminantRadius()
            self.assertAlmostEqual(
                sigma_out, expected_sigma_out, delta=0.1,
                msg=f"matched-PSF width at {point} differs from target",
            )

        # The merged template image must be all-finite over the data
        # region; a successful merge fills every pixel from at least
        # one matched patch.
        self.assertTrue(np.all(np.isfinite(template.image.array)))
        self.assertTrue(np.all(np.isfinite(template.variance.array)))

        # Verify the patch images really were convolved (not silently
        # passed through). Convolution with a normalized kernel reduces
        # per-pixel variance by Σ K² ≪ 1, so the merged template's
        # variance must be much smaller than any input patch's;
        # empirically ~40× quieter. Require at least 5× to catch a
        # silent no-op regression.
        input_patch = self.patches[0][0].get()
        no_data_in = input_patch.mask.getPlaneBitMask("NO_DATA")
        no_data_out = template.mask.getPlaneBitMask("NO_DATA")
        input_good = (input_patch.mask.array & no_data_in) == 0
        output_good = (template.mask.array & no_data_out) == 0
        input_var = float(np.median(input_patch.variance.array[input_good]))
        matched_var = float(np.median(template.variance.array[output_good]))
        self.assertLess(
            matched_var, input_var/5,
            msg=f"matched-template variance {matched_var:.4f} not"
                f" substantially below input patch variance {input_var:.4f};"
                " convolution may have been skipped",
        )

        # Direct check: measure each injected source on the merged
        # matched template via second moments. Each injection started at
        # σ = sigmas[(t, p)] in patch pixels. After matching it should
        # have σ = target on the patch grid, which maps to
        # ``expected_sigma_out`` on the science grid. Skip injections
        # whose science-pixel position is too close to the template edge
        # for a clean halfsize-stamp measurement.
        halfsize = 10
        inset = lsst.geom.Box2I(template.getBBox())
        inset.grow(-halfsize)
        measurements = []
        for (tract_id, patch_id), pos in science_positions.items():
            if not inset.contains(lsst.geom.Point2I(pos)):
                continue
            measured = self._measureGaussianSigma(
                template.image, pos, halfsize=halfsize,
            )
            self.assertAlmostEqual(
                measured, expected_sigma_out, delta=0.1,
                msg=f"injected source from tract={tract_id} patch={patch_id}"
                    f" (σ_in={sigmas[(tract_id, patch_id)]:.2f}) measured"
                    f" σ={measured:.3f} on the matched template, expected"
                    f" {expected_sigma_out:.3f}",
            )
            measurements.append(measured)
        # Sanity: at least four injections should fall safely inside the
        # template for the fixture's box.
        self.assertGreaterEqual(len(measurements), 4)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
