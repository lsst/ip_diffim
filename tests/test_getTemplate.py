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
from unittest import mock
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

import lsst.geom as geom
import lsst.afw.image
import lsst.afw.math
from lsst.afw.image import ExposureF, MaskedImageF
from lsst.daf.butler import DataCoordinate, DimensionUniverse
import lsst.ip.diffim
from lsst.ip.diffim.getTemplate import (
    GetTemplateTask,
    GetDcrTemplateTask,
    GetTemplateConnections,
    GetDcrTemplateConnections,
)
import lsst.meas.algorithms
import lsst.meas.base.tests
import lsst.pipe.base as pipeBase
from lsst.pipe.base import NoWorkFound
import lsst.skymap
import lsst.utils.tests

from contextlib import contextmanager

from utils import generate_data_id

# Change this to True, `setup display_ds9`, and open ds9 (or use another afw
# display backend) to show the tract/patch layouts on the image.
debug = False
if debug:
    import lsst.afw.display
    display = lsst.afw.display.Display()
    display.frame = 1


@contextmanager
def mock_polygon_intersection(intersects=True, area=100):
    """Context manager to patch lsst.afw.geom.Polygon with controlled
    intersection behavior.

    Parameters
    ----------
    intersects : bool
        If True, polygons intersect; if False, no intersection.
    area : float
        Area returned by `calculateArea` for intersectionSingle.
    """
    with patch("lsst.afw.geom.Polygon") as mockPolygon:
        poly_instance = MagicMock()
        if intersects:
            poly_instance.intersection.return_value = True
            poly_instance.intersectionSingle.return_value.calculateArea.return_value = area
        else:
            poly_instance.intersection.return_value = None
            poly_instance.intersectionSingle.return_value.calculateArea.return_value = 0
        mockPolygon.return_value = poly_instance
        yield


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
        self.getTemplateTask = lsst.ip.diffim.getTemplate.GetTemplateTask
        self.expectedNoDataSum = 20990
        self.scale = 0.2  # arcsec/pixel
        self.skymap = self._makeSkymap()
        self.patches = collections.defaultdict(list)
        self.dataIds = collections.defaultdict(list)
        self.exposure = self._makeExposure()
        self.useDcr = False

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
                    tract=tract.tract_id,
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
        # self.useDcr set in the right places, True & False in setup
        if self.useDcr:
            self.assertEqual(template.psf.getComponentCount(), nPsfs*self.dcrNumSubfilters)
        else:
            self.assertEqual(template.psf.getComponentCount(), nPsfs)

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
        # Variance plane ==2 in the original image, but the warped images will
        # have some structure due to the warping.
        self.assertImagesAlmostEqual(template.variance, self.exposure[expectedBox].variance,
                                     rtol=0.55, msg="variance planes differ")
        # Not checking the mask, as warping changes the sizes of the masks.

    def testRunQuantum(self):
        """
        Test runQuantum (non-DCR version).
        """
        # Create a local task instance
        config = GetTemplateTask.ConfigClass()
        task = GetTemplateTask(config=config)

        # Create mock inputs
        mockCoaddExposures = ["coadd1", "coadd2"]
        mockWcs = mock.Mock(name="wcs")
        mockBBox = mock.Mock(name="bbox")
        mockSkymap = mock.Mock(name="skymap")

        # Bundle into dictionary that runQuantum expects
        inputDict = {
            "coaddExposures": mockCoaddExposures,
            "wcs": mockWcs,
            "bbox": mockBBox,
            "skyMap": mockSkymap,
        }

        # Mock Butler QC object
        butlerQC = mock.Mock()
        butlerQC.get.return_value = inputDict
        butlerQC.quantum.dataId = {"physical_filter": "r"}
        butlerQC.put = mock.Mock()

        # Mock getExposures to return a simple object with expected attributes
        mockResults = mock.Mock()
        mockResults.coaddExposures = ["coaddExposure1", "coaddExposure2"]
        mockResults.dataIds = ["dataId1", "dataId2"]
        task.getExposures = mock.Mock(return_value=mockResults)

        # Mock the run method to just return a dummy output
        task.run = mock.Mock(return_value="final_output")

        # Call runQuantum
        outputRefs = "outputRefs"
        inputRefs = "inputRefs"
        task.runQuantum(butlerQC, inputRefs, outputRefs)

        # Check that get was called correctly
        butlerQC.get.assert_called_once_with(inputRefs)

        # Check that getExposures was called with correct arguments
        task.getExposures.assert_called_once_with(
            mockCoaddExposures, mockBBox, mockSkymap, mockWcs
        )

        # Check that run was called with the right parameters
        task.run.assert_called_once_with(
            coaddExposureHandles=mockResults.coaddExposures,
            bbox=mockBBox,
            wcs=mockWcs,
            dataIds=mockResults.dataIds,
            physical_filter="r",
        )

        # Check that put was called with the outputs
        butlerQC.put.assert_called_once_with("final_output", outputRefs)

    def testGetExposuresRaisesWhenWcsIsNone(self):
        """
        Test getExposures raises NoWorkFound when WCS is None.
        """
        # Create task
        config = GetTemplateTask.ConfigClass()
        task = GetTemplateTask(config=config)

        # Mock inputs
        coaddHandles = [MagicMock()]
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Point2I(100, 100))
        wcs = None  # deliberately None

        with self.assertRaises(NoWorkFound):
            task.getExposures(coaddHandles, bbox, self.skymap, wcs)

    def testGetExposuresRaisesWhenNoOverlap(self):
        """
        Test getExposures raises NoWorkFound when no patches overlap the
        detector bbox.
        """
        task = GetTemplateTask(config=GetTemplateTask.ConfigClass())
        bbox = geom.Box2I(geom.Point2I(100, 100), geom.Point2I(200, 200))

        # Mock patch, tract, skymap
        mock_patch = MagicMock()
        mock_patch.getOuterBBox.return_value = geom.Box2I(geom.Point2I(0, 0), geom.Point2I(10, 10))
        mock_tract = MagicMock()
        mock_tract.getWcs.return_value = MagicMock()
        mock_tract.__getitem__.return_value = mock_patch
        skymap = {0: mock_tract}

        coadd_handle = MagicMock()
        coadd_handle.dataId = {"tract": 0, "patch": 0, "subfilter": 0}
        coaddHandles = [coadd_handle]

        with mock_polygon_intersection(intersects=False):
            with self.assertRaises(NoWorkFound):
                task.getExposures(coaddHandles, bbox, skymap, MagicMock())

    def testGetExposuresWithOverlap(self):
        """
        Test getExposures with an overlapping patch exercises the intersection
        branch.
        """
        task = GetTemplateTask(config=GetTemplateTask.ConfigClass())
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Point2I(100, 100))

        # Mock patch and tract
        overlap_patch = MagicMock()
        overlap_patch.getOuterBBox.return_value = geom.Box2I(geom.Point2I(10, 10), geom.Point2I(50, 50))
        mock_tract = MagicMock()
        mock_tract.getWcs.return_value = MagicMock()
        mock_tract.__getitem__.return_value = overlap_patch
        skymap = {0: mock_tract}

        coadd_handle = MagicMock()
        coadd_handle.dataId = {"tract": 0, "patch": 0, "subfilter": 0}
        coaddHandles = [coadd_handle]

        # Use the mock_polygon_intersection helper
        with mock_polygon_intersection(intersects=True, area=100):
            result = task.getExposures(coaddHandles, bbox, skymap, MagicMock())

        # Check that the overlapping exposure is returned
        self.assertIn(0, result.coaddExposures)
        self.assertIn(0, result.dataIds)
        self.assertEqual(result.coaddExposures[0][0], coadd_handle)
        self.assertEqual(result.dataIds[0][0], coadd_handle.dataId)

    def testCheckHighVarianceLowGoodFraction(self):
        """
        Test that checkHighVariance logs a message when
        goodFraction < highVarianceMaskFraction.
        """
        # Create a task with a known threshold
        config = GetTemplateTask.ConfigClass()
        config.highVarianceMaskFraction = 0.8  # Require 80% good pixels
        task = GetTemplateTask(config=config)

        # Mock the template exposure
        mock_mask = mock.Mock()
        # Simulate a small mask array where only 50% are "good"
        mask_array = np.array([[0, 1], [0, 1]], dtype=int)
        mock_mask.array = mask_array

        # The ignored mask bits — pretend bit 1 is ignored
        mock_mask.getPlaneBitMask.return_value = 1

        # addMaskPlane returns a bit number for "HIGH_VARIANCE"
        mock_mask.addMaskPlane.return_value = 32

        # Attach the mask to a fake exposure
        template = mock.Mock()
        template.mask = mock_mask

        # Mock the varianceBackground.config.ignoredPixelMask used internally
        task.varianceBackground = mock.Mock()
        task.varianceBackground.config.ignoredPixelMask = ["BAD", "NO_DATA"]

        # Patch the task logger to capture info messages
        with mock.patch.object(task.log, "info") as mock_log:
            task.checkHighVariance(template)

        # Assert that warning is logged when goodFraction = 0.5 (50%) and
        # threshold = 0.8 (80%)
        mock_log.assert_called_once()
        log_args = mock_log.call_args[0][0]
        assert "Not setting HIGH_VARIANCE mask plane" in log_args

        # addMaskPlane should still be called
        mock_mask.addMaskPlane.assert_called_once_with("HIGH_VARIANCE")

        # getPlaneBitMask should be called with the ignored mask list
        mock_mask.getPlaneBitMask.assert_called_once_with(
            task.varianceBackground.config.ignoredPixelMask
        )

    def testCheckInputsMultipleBands(self):
        """
        Test that _checkInputs raises RuntimeError when multiple bands are
        found in dataIds.
        """
        task = GetTemplateTask(config=GetTemplateTask.ConfigClass())

        # Two tracts with conflicting band entries
        dataIds = {
            1: [{"band": "g"}],
            2: [{"band": "r"}],
        }

        # Instantiate coadd exposures
        coaddExposures = {1: [], 2: []}

        # Check that the RuntimeError is raised
        with pytest.raises(RuntimeError, match="multiple bands"):
            task._checkInputs(dataIds, coaddExposures)

    def testCheckInputsDifferentPhotoCalibs(self):
        """
        Test that _checkInputs raises a RuntimeError when coadd exposures have
        different photoCalibs.
        """

        task = GetTemplateTask(config=GetTemplateTask.ConfigClass())

        # All dataIds have same band
        dataIds = {1: [{"band": "r"}]}

        # Create mock exposures that return different photoCalibs
        exp1 = mock.Mock()
        exp2 = mock.Mock()
        exp1.get.return_value = "photoCalib1"
        exp2.get.return_value = "photoCalib2"

        coaddExposures = {1: [exp1, exp2]}

        # Check that the correct RuntimeError is raised
        with pytest.raises(RuntimeError, match="different photoCalibs"):
            task._checkInputs(dataIds, coaddExposures)

    def testRunOneTractInput(self):
        """Test a bounding box that fully fits inside one tract, with only
        that tract passed as input. This checks that the code handles a single
        tract input correctly.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = self.getTemplateTask()
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
        task = self.getTemplateTask()
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
        task = self.getTemplateTask()
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
        task = self.getTemplateTask()
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
        task = self.getTemplateTask()
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        no_data = (result.template.mask.array & result.template.mask.getPlaneBitMask("NO_DATA")) != 0
        self.assertTrue(np.isfinite(result.template.image.array).all())
        self.assertTrue(np.isfinite(result.template.variance.array).all())
        if self.useDcr:
            self.assertEqual(no_data.sum(), 0)
        else:
            self.assertEqual(no_data.sum(), self.expectedNoDataSum)

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
        task = self.getTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        if debug:
            _showTemplate(box, result.template)
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 9)
        # We just check that the pixel values are all finite. We cannot check that pixel values
        # in the template are closer to the original anymore.
        self.assertTrue(np.isfinite(result.template.image.array).all())


class GetDcrTemplateTaskTestCase(GetTemplateTaskTestCase):
    """Test that GetDcrTemplateTask works on both one tract and multiple tract
    input coadd exposures.

    Makes a synthetic exposure large enough to fit four small tracts with 2x2
    (300x300 pixel) patches each, extracts pixels for those patches by warping,
    and tests GetDCRTemplateTask's output against boxes that overlap various
    combinations of one or multiple tracts.
    """

    def setUp(self):
        self.dcrNumSubfilters = 3
        super().setUp()
        getTemplateTask = lsst.ip.diffim.getTemplate.GetDcrTemplateTask
        config = getTemplateTask.ConfigClass()
        config.bandwidth = 147.0
        config.effectiveWavelength = 478.5

        self.useDcr = True
        self.expectedNoDataSum = 0
        # self.bandwidth = 147.0
        config.effectiveWavelength = 478.5
        # self.effectiveWavelength = 478.5
        getTemplateTask = getTemplateTask(config=config)
        # self.getTemplateTask = getTemplateTask(config=config)

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

            # This is mostly stolen from pipe_tasks warpAndPsfMatch, but
            # ip_diffim cannot depend on pipe_tasks.
            xyTransform = lsst.afw.geom.makeWcsPairTransform(self.exposure.wcs, patch.wcs)
            warpedPsf = lsst.meas.algorithms.WarpedPsf(self.exposure.psf, xyTransform)
            warped = warper.warpExposure(patch.wcs, self.exposure, destBBox=box)
            warped.setPsf(warpedPsf)
            for subfilter in range(self.dcrNumSubfilters):
                dataId = generate_data_id(
                    tract=tract.tract_id,
                    patch=patch.getSequentialIndex(),
                    subfilter=subfilter
                )
                dataRef = pipeBase.InMemoryDatasetHandle(
                    warped,
                    storageClass="ExposureF",
                    copy=True,
                    dataId=dataId
                )
                self.patches[tract.tract_id].append(dataRef)
                dataCoordinate = DataCoordinate.standardize({"tract": tract.tract_id,
                                                             "patch": patch.getSequentialIndex(),
                                                             "band": "a",
                                                             "subfilter": subfilter,
                                                             "skymap": "skymap"},
                                                            universe=DimensionUniverse())
                self.dataIds[tract.tract_id].append(dataCoordinate)

    def testConnectionsRemovesCoaddExposures(self):
        """
        Test that GetDcrTemplateConnections.__init__ removes 'coaddExposures'
        from the list of inputs inherited from GetTemplateConnections.
        """
        # Instantiate both the base and subclass with their configs
        base_config = GetTemplateTask.ConfigClass()
        base_connections = GetTemplateConnections(config=base_config)

        dcr_config = GetDcrTemplateTask.ConfigClass()
        dcr_connections = GetDcrTemplateConnections(config=dcr_config)

        # Sanity check: base class *should* have "coaddExposures"
        assert "coaddExposures" in base_connections.inputs, (
            "Sanity check failed: base GetTemplateConnections should include 'coaddExposures'"
        )

        # Verify that the subclass removes it
        assert "coaddExposures" not in dcr_connections.inputs, (
            "GetDcrTemplateConnections.__init__ should remove 'coaddExposures' from inputs"
        )

        # Optional: check that DCR-specific inputs exist
        for expected in ("dcrCoadds", "visitInfo"):
            assert expected in dcr_connections.inputs, f"Expected '{expected}' missing in DCR connections"

    def testDcrRunQuantum(self):
        """
        Test dcr version of runQuantum.
        """

        # Create a local task instance
        config = GetDcrTemplateTask.ConfigClass()
        config.bandwidth = 147.0
        config.effectiveWavelength = 478.5
        task = GetDcrTemplateTask(config=config)

        # Create mock inputs
        mockDcrCoadds = ["coadd1", "coadd2"]
        mockWcs = mock.Mock(name="wcs")
        mockBBox = mock.Mock(name="bbox")
        mockSkymap = mock.Mock(name="skymap")
        mockVisitInfo = mock.Mock(name="visitInfo")

        # Bundle into the dictionary that runQuantum expects
        inputDict = {
            "dcrCoadds": mockDcrCoadds,
            "wcs": mockWcs,
            "bbox": mockBBox,
            "skyMap": mockSkymap,
            "visitInfo": mockVisitInfo,
        }

        # Mock Butler QC object
        butlerQC = mock.Mock()
        butlerQC.get.return_value = inputDict
        butlerQC.quantum.dataId = {"physical_filter": "r"}
        butlerQC.put = mock.Mock()

        # Mock getExposures to return a simple object with expected attributes
        mockResults = mock.Mock()
        mockResults.coaddExposures = ["coaddExposure1", "coaddExposure2"]
        mockResults.dataIds = ["dataId1", "dataId2"]
        task.getExposures = mock.Mock(return_value=mockResults)

        # Mock the run method to just return a dummy output
        task.run = mock.Mock(return_value="final_output")

        # Call runQuantum
        outputRefs = "outputRefs"
        inputRefs = "inputRefs"
        task.runQuantum(butlerQC, inputRefs, outputRefs)

        # Check that get was called correctly
        butlerQC.get.assert_called_once_with(inputRefs)

        # Check that getExposures was called with correct arguments
        task.getExposures.assert_called_once_with(
            mockDcrCoadds, mockBBox, mockSkymap, mockWcs, mockVisitInfo
        )

        # Check that run was called with the right parameters
        task.run.assert_called_once_with(
            coaddExposureHandles=mockResults.coaddExposures,
            bbox=mockBBox,
            wcs=mockWcs,
            dataIds=mockResults.dataIds,
            physical_filter="r",
        )

        # Check that put was called with the outputs
        butlerQC.put.assert_called_once_with("final_output", outputRefs)

    def testValidate(self):
        # Check that function does not raise ValueError when
        # self.effectiveWavelength and self.bandwidth are not None.
        self.getTemplateTask()

        # Check that function raises ValueError if self.effectiveWavelength is
        # None or self.bandwidth is None
        with self.assertRaises(ValueError):
            task = lsst.ip.diffim.getTemplate.GetDcrTemplateTask
            valueConfig = task.ConfigClass()
            task(config=valueConfig)

    def testDcrGetExposuresRaisesWhenWcsIsNone(self):
        """
        Test that getExposures raises NoWorkFound when WCS is None.
        """
        # Create task with valid DCR config
        config = GetDcrTemplateTask.ConfigClass()
        config.bandwidth = 147.0
        config.effectiveWavelength = 478.5
        config.numSubfilters = 3
        task = GetDcrTemplateTask(config=config)

        # Mock inputs
        coaddHandles = [mock.Mock(name="handle")]
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Point2I(100, 100))
        skymap = mock.Mock(name="skymap")
        wcs = None  # deliberately None
        visitInfo = mock.Mock(name="visitInfo")  # required argument

        # Expect NoWorkFound when wcs is None
        with self.assertRaises(NoWorkFound):
            task.getExposures(coaddHandles, bbox, skymap, wcs, visitInfo)

    def testDcrGetExposuresRaisesWhenNoOverlap(self):
        """
        Test that getExposures raises NoWorkFound when no patches overlap the
        detector bbox.
        """
        # Create task with valid DCR config
        config = GetDcrTemplateTask.ConfigClass()
        config.bandwidth = 147.0
        config.effectiveWavelength = 478.5
        config.numSubfilters = 3
        task = GetDcrTemplateTask(config=config)

        # Mock coadd exposure handles and their dataIds
        coaddRef = mock.MagicMock()
        coaddRef.dataId = {"tract": 123, "patch": "1,1", "subfilter": 0}
        coaddHandles = [coaddRef]

        # Mock patch
        mock_patch = mock.MagicMock()
        mock_patch.getOuterBBox.return_value = geom.Box2I(
            geom.Point2I(0, 0), geom.Point2I(10, 10)
        )

        # Mock tract
        mock_tract = mock.MagicMock()
        mock_tract.getWcs.return_value = mock.MagicMock()
        mock_tract.__getitem__.return_value = mock_patch

        # Mock skymap
        skymap = {123: mock_tract}

        # Mock WCS that returns non-overlapping polygons
        mock_wcs = mock.MagicMock()
        # Polygon intersection will be empty → forces no overlap
        with mock_polygon_intersection(intersects=False):
            # Non-overlapping detector bounding box
            bbox = geom.Box2I(geom.Point2I(100, 100), geom.Point2I(200, 200))

            # Mock visitInfo to satisfy the argument
            visitInfo = mock.MagicMock()

            with self.assertRaises(NoWorkFound):
                task.getExposures(coaddHandles, bbox, skymap, mock_wcs, visitInfo)

    def testDcrGetExposuresDataValidity(self):
        """Test getExposures includes only overlapping patches and that image,
        mask, and variance are valid."""
        # Configure DCR task
        config = GetDcrTemplateTask.ConfigClass()
        config.bandwidth = 147.0  # Required for DCR
        config.effectiveWavelength = 478.5  # Required for DCR
        config.numSubfilters = 3  # Number of DCR subfilters expected per patch
        task = GetDcrTemplateTask(config=config)

        # Define detector bounding box
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Point2I(100, 100))

        # Mock WCS
        mockWcs = MagicMock()
        mockWcs.skyToPixel.side_effect = lambda corners: corners

        def make_mock_exposure(shape=(50, 50), fill_value=100):
            mi = MaskedImageF(shape[0], shape[1])
            mi.image.array[:] = fill_value
            mi.variance.array[:] = 1.0
            mi.mask.array[:] = 0
            return ExposureF(mi)

        # Mock patches
        # One patch overlapping detector
        overlap_patch = MagicMock()
        overlap_patch.getOuterBBox.return_value = geom.Box2I(geom.Point2I(10, 10),
                                                             geom.Point2I(50, 50))
        # One patch outside detector
        non_overlap_patch = MagicMock()
        non_overlap_patch.getOuterBBox.return_value = geom.Box2I(geom.Point2I(200, 200),
                                                                 geom.Point2I(250, 250))

        # Mock tracts containing patches
        overlap_tract = MagicMock()
        overlap_tract.getWcs.return_value = mockWcs
        overlap_tract.__getitem__.side_effect = lambda patch_id: overlap_patch

        non_overlap_tract = MagicMock()
        non_overlap_tract.getWcs.return_value = mockWcs
        non_overlap_tract.__getitem__.side_effect = lambda patch_id: non_overlap_patch

        # Mock SkyMap
        skymap = MagicMock()
        skymap.__getitem__.side_effect = lambda tract_id: overlap_tract if tract_id == 0 else \
            non_overlap_tract

        # Create mock exposure handles for all subfilters
        exposureHandles = []
        for subfilter in range(config.numSubfilters):
            overlap_handle = MagicMock()
            overlap_handle.dataId = {"tract": 0, "patch": 0, "subfilter": subfilter}
            non_overlap_handle = MagicMock()
            non_overlap_handle.dataId = {"tract": 1, "patch": 0, "subfilter": subfilter}
            exposureHandles.extend([overlap_handle, non_overlap_handle])

        # Patch getDcrModel to return a real masked image for tract 0
        mock_exposure = make_mock_exposure()
        task.getDcrModel = MagicMock(return_value={0: mock_exposure})

        # Patch afwGeom.Polygon to control intersection behavior
        with patch("lsst.afw.geom.Polygon") as mockPolygon:
            poly_instance = MagicMock()

            # Only tract 0 patches intersect the detector
            def intersection(other):
                if other == geom.Box2D(bbox):
                    return True  # intersection exists
                return None

            poly_instance.intersection.side_effect = intersection
            poly_instance.intersectionSingle.return_value.calculateArea.return_value = 100
            mockPolygon.return_value = poly_instance

            # Call getExposures
            result = task.getExposures(exposureHandles, bbox, skymap, mockWcs, visitInfo=MagicMock())

        # Validate which tracts were included
        included_tracts = set(result.coaddExposures.keys())
        self.assertIn(0, included_tracts, "Overlapping tract should be included")
        self.assertNotIn(1, included_tracts, "Non-overlapping tract should not be included")

        # Validate image, variance, and mask data
        exp = result.coaddExposures[0]  # only overlapping tract
        mi = exp.maskedImage

        # Check image has non-zero values in overlap region
        self.assertTrue(np.any(mi.image.array != 0), "No image data in overlapping region")

        # Check variance is set
        self.assertTrue(np.any(mi.variance.array > 0), "Variance not set in overlapping region")

        # Check mask plane is clean (no unexpected bits set)
        self.assertTrue(np.all(mi.mask.array == 0), "Mask plane has unexpected bits set")

    def testGetDcrModelReturnsExposures(self):
        """
        Test that getDcrModel returns a dict of exposures per tract, using
        mocked DcrModel.
        """
        # Create a task with minimal valid config
        config = GetDcrTemplateTask.ConfigClass()
        config.bandwidth = 150.0
        config.effectiveWavelength = 500.0
        config.numSubfilters = 3
        task = GetDcrTemplateTask(config=config)

        visitInfo = MagicMock()
        coaddRefs = [MagicMock(), MagicMock()]
        patchList = {0: [0]}

        with patch("lsst.ip.diffim.getTemplate._selectDataRef", return_value=True):
            # Create a mock instance to be returned by fromQuantum
            mock_dcrModel_instance = MagicMock()
            mock_dcrModel_instance.buildMatchedExposureHandle.return_value = "mock_exposure"

            # Patch the classmethod fromQuantum to return mock instance
            with patch("lsst.ip.diffim.getTemplate.DcrModel.fromQuantum",
                       return_value=mock_dcrModel_instance) as mock_fromQuantum:
                result = task.getDcrModel(patchList, coaddRefs, visitInfo)

        # Verify returned dict
        assert isinstance(result, dict)
        assert 0 in result
        assert result[0] == ["mock_exposure"]

        # Check buildMatchedExposureHandle called with visitInfo
        mock_dcrModel_instance.buildMatchedExposureHandle.assert_called_with(visitInfo=visitInfo)

        # Check fromQuantum was called with expected arguments
        mock_fromQuantum.assert_called_with(
            coaddRefs,
            config.effectiveWavelength,
            config.bandwidth,
            config.numSubfilters
        )

    def testCheckPatchList(self):
        """
        Test that checkPatchList correctly validates the number of DCR
        subfilters per patch.

        Setup:
        - For each tract in self.patches, create a list of patch handles.
        - Each handle is repeated once for each expected DCR subfilter.
          This simulates having multiple DCR-matched exposures per patch,
          one for each subfilter.

        Expected behavior:
        - If the number of patch entries matches task.config.numSubfilters,
          checkPatchList passes.
        - If a patch has too few or too many entries, checkPatchList raises a
          RuntimeError.
        """

        # Create a temporary instance of the task
        taskClass = lsst.ip.diffim.getTemplate.GetDcrTemplateTask
        config = taskClass.ConfigClass()
        config.bandwidth = 147.0
        config.effectiveWavelength = 478.5
        config.numSubfilters = 3
        task = taskClass(config=config)

        # 1. Prepare a patch list with correct number of subfilters
        patchList = {}
        numSubfilters = task.config.numSubfilters
        for tract, handles in self.patches.items():
            # Repeat each handle once per subfilter
            patchList[tract] = [handle for handle in handles for _ in range(numSubfilters)]

        # Should not raise an exception
        task.checkPatchList(patchList)

        # 2. Test failure case: remove one entry for a patch
        badPatchList = {tract: ids[:-1] for tract, ids in patchList.items()}
        with self.assertRaises(RuntimeError):
            task.checkPatchList(badPatchList)

    def testSelectDataRef(self):
        """Test that the input tract and patch match the those from the dataId.
        """
        coaddRef1 = self.patches[0][0]
        coaddRef2 = self.patches[0][4]
        tract = coaddRef1.dataId['tract']
        patch = coaddRef1.dataId['patch']
        selectDataRef = lsst.ip.diffim.getTemplate._selectDataRef(coaddRef1, tract, patch)
        self.assertTrue(selectDataRef)

        patch2 = coaddRef2.dataId['patch']
        selectDataRef2 = lsst.ip.diffim.getTemplate._selectDataRef(coaddRef1, tract, patch2)
        self.assertFalse(selectDataRef2)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
