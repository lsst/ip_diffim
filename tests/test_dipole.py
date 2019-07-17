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

import unittest

import numpy as np

import lsst.utils.tests
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.ip.diffim as ipDiffim

display = False
try:
    display
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)

sigma2fwhm = 2.*np.sqrt(2.*np.log(2.))


def makePluginAndCat(alg, name, control, metadata=False, centroid=None):
    schema = afwTable.SourceTable.makeMinimalSchema()
    if centroid:
        schema.addField(centroid + "_x", type=float)
        schema.addField(centroid + "_y", type=float)
        schema.addField(centroid + "_flag", type='Flag')
        schema.getAliasMap().set("slot_Centroid", centroid)
    if metadata:
        plugin = alg(control, name, schema, dafBase.PropertySet())
    else:
        plugin = alg(control, name, schema)
    cat = afwTable.SourceCatalog(schema)
    return plugin, cat


def createDipole(w, h, xc, yc, scaling=100.0, fracOffset=1.2):
    # Make random noise image: set image plane to normal distribution
    image = afwImage.MaskedImageF(w, h)
    image.set(0)
    array = image.getImage().getArray()
    array[:, :] = np.random.randn(w, h)
    # Set variance to 1.0
    var = image.getVariance()
    var.set(1.0)

    if display:
        afwDisplay.Display(frame=1).mtv(image, title="Original image")
        afwDisplay.Display(frame=2).mtv(image.getVariance(), title="Original variance")

    # Create Psf for dipole creation and measurement
    psfSize = 17
    psf = measAlg.DoubleGaussianPsf(psfSize, psfSize, 2.0, 3.5, 0.1)
    psfFwhmPix = sigma2fwhm*psf.computeShape().getDeterminantRadius()
    psfim = psf.computeImage().convertF()
    psfim *= scaling/psf.computePeak()
    psfw, psfh = psfim.getDimensions()
    psfSum = np.sum(psfim.getArray())

    # Create the dipole, offset by fracOffset of the Psf FWHM (pixels)
    offset = fracOffset*psfFwhmPix//2
    array = image.getImage().getArray()
    xp = int(xc - psfw//2 + offset)
    yp = int(yc - psfh//2 + offset)
    array[yp:yp + psfh, xp:xp + psfw] += psfim.getArray()

    xn = int(xc - psfw//2 - offset)
    yn = int(yc - psfh//2 - offset)
    array[yn:yn + psfh, xn:xn + psfw] -= psfim.getArray()

    if display:
        afwDisplay.Display(frame=3).mtv(image, title="With dipole")

    # Create an exposure, detect positive and negative peaks separately
    exp = afwImage.makeExposure(image)
    exp.setPsf(psf)
    config = measAlg.SourceDetectionConfig()
    config.thresholdPolarity = "both"
    config.reEstimateBackground = False
    schema = afwTable.SourceTable.makeMinimalSchema()
    task = measAlg.SourceDetectionTask(schema, config=config)
    table = afwTable.SourceTable.make(schema)
    results = task.makeSourceCatalog(table, exp)
    if display:
        afwDisplay.Display(frame=4).mtv(image, title="Detection plane")

    # Merge them together
    assert(len(results.sources) == 2)
    fpSet = results.fpSets.positive
    fpSet.merge(results.fpSets.negative, 0, 0, False)
    sources = afwTable.SourceCatalog(table)
    fpSet.makeSources(sources)
    assert(len(sources) == 1)
    s = sources[0]
    assert(len(s.getFootprint().getPeaks()) == 2)

    return psf, psfSum, exp, s


class DipoleAlgorithmTest(lsst.utils.tests.TestCase):
    """ A test case for dipole algorithms"""

    def setUp(self):
        np.random.seed(666)
        self.w, self.h = 100, 100  # size of image
        self.xc, self.yc = 50, 50  # location of center of dipole

    def testNaiveDipoleCentroid(self):
        control = ipDiffim.DipoleCentroidControl()
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        plugin, cat = makePluginAndCat(ipDiffim.NaiveDipoleCentroid, "test", control, centroid="centroid")
        source = cat.addNew()
        source.set("centroid_x", 50)
        source.set("centroid_y", 50)
        source.setFootprint(s.getFootprint())
        plugin.measure(source, exposure)
        for key in ("_pos_x", "_pos_y", "_pos_xErr", "_pos_yErr", "_pos_flag",
                    "_neg_x", "_neg_y", "_neg_xErr", "_neg_yErr", "_neg_flag"):
            try:
                source.get("test" + key)
            except Exception:
                self.fail()

    def testNaiveDipoleFluxControl(self):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        control = ipDiffim.DipoleFluxControl()
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        plugin, cat = makePluginAndCat(ipDiffim.NaiveDipoleFlux, "test", control, centroid="centroid")
        source = cat.addNew()
        source.set("centroid_x", 50)
        source.set("centroid_y", 50)
        source.setFootprint(s.getFootprint())
        plugin.measure(source, exposure)
        for key in ("_pos_instFlux", "_pos_instFluxErr", "_pos_flag", "_npos",
                    "_neg_instFlux", "_neg_instFluxErr", "_neg_flag", "_nneg"):
            try:
                source.get("test" + key)
            except Exception:
                self.fail()

    def testPsfDipoleFluxControl(self):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        control = ipDiffim.PsfDipoleFluxControl()
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        plugin, cat = makePluginAndCat(ipDiffim.PsfDipoleFlux, "test", control, centroid="centroid")
        source = cat.addNew()
        source.set("centroid_x", 50)
        source.set("centroid_y", 50)
        source.setFootprint(s.getFootprint())
        plugin.measure(source, exposure)
        for key in ("_pos_instFlux", "_pos_instFluxErr", "_pos_flag",
                    "_neg_instFlux", "_neg_instFluxErr", "_neg_flag"):
            try:
                source.get("test" + key)
            except Exception:
                self.fail()

    def testAll(self):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        self.measureDipole(s, exposure)

    def _makeModel(self, exposure, psf, fp, negCenter, posCenter):

        negPsf = psf.computeImage(negCenter).convertF()
        posPsf = psf.computeImage(posCenter).convertF()
        negPeak = psf.computePeak(negCenter)
        posPeak = psf.computePeak(posCenter)
        negPsf /= negPeak
        posPsf /= posPeak

        model = afwImage.ImageF(fp.getBBox())
        negModel = afwImage.ImageF(fp.getBBox())
        posModel = afwImage.ImageF(fp.getBBox())

        # The center of the Psf should be at negCenter, posCenter
        negPsfBBox = negPsf.getBBox()
        posPsfBBox = posPsf.getBBox()
        modelBBox = model.getBBox()

        # Portion of the negative Psf that overlaps the montage
        negOverlapBBox = geom.Box2I(negPsfBBox)
        negOverlapBBox.clip(modelBBox)
        self.assertFalse(negOverlapBBox.isEmpty())

        # Portion of the positivePsf that overlaps the montage
        posOverlapBBox = geom.Box2I(posPsfBBox)
        posOverlapBBox.clip(modelBBox)
        self.assertFalse(posOverlapBBox.isEmpty())

        negPsfSubim = type(negPsf)(negPsf, negOverlapBBox)
        modelSubim = type(model)(model, negOverlapBBox)
        negModelSubim = type(negModel)(negModel, negOverlapBBox)
        modelSubim += negPsfSubim  # just for debugging
        negModelSubim += negPsfSubim  # for fitting

        posPsfSubim = type(posPsf)(posPsf, posOverlapBBox)
        modelSubim = type(model)(model, posOverlapBBox)
        posModelSubim = type(posModel)(posModel, posOverlapBBox)
        modelSubim += posPsfSubim
        posModelSubim += posPsfSubim

        data = afwImage.ImageF(exposure.getMaskedImage().getImage(), fp.getBBox())
        var = afwImage.ImageF(exposure.getMaskedImage().getVariance(), fp.getBBox())
        matrixNorm = 1./np.sqrt(np.median(var.getArray()))

        if display:
            afwDisplay.Display(frame=5).mtv(model, title="Unfitted model")
            afwDisplay.Display(frame=6).mtv(data, title="Data")

        posPsfSum = np.sum(posPsf.getArray())
        negPsfSum = np.sum(negPsf.getArray())

        M = np.array((np.ravel(negModel.getArray()), np.ravel(posModel.getArray()))).T.astype(np.float64)
        B = np.array((np.ravel(data.getArray()))).astype(np.float64)
        M *= matrixNorm
        B *= matrixNorm

        # Numpy solution
        fneg0, fpos0 = np.linalg.lstsq(M, B, rcond=-1)[0]

        # Afw solution
        lsq = afwMath.LeastSquares.fromDesignMatrix(M, B, afwMath.LeastSquares.DIRECT_SVD)
        fneg, fpos = lsq.getSolution()

        # Should be exaxtly the same as each other
        self.assertAlmostEqual(1e-2*fneg0, 1e-2*fneg)
        self.assertAlmostEqual(1e-2*fpos0, 1e-2*fpos)

        # Recreate model
        fitted = afwImage.ImageF(fp.getBBox())
        negFit = type(negPsf)(negPsf, negOverlapBBox, afwImage.PARENT, True)
        negFit *= float(fneg)
        posFit = type(posPsf)(posPsf, posOverlapBBox, afwImage.PARENT, True)
        posFit *= float(fpos)

        fitSubim = type(fitted)(fitted, negOverlapBBox)
        fitSubim += negFit
        fitSubim = type(fitted)(fitted, posOverlapBBox)
        fitSubim += posFit
        if display:
            afwDisplay.Display(frame=7).mtv(fitted, title="Fitted model")

        fitted -= data

        if display:
            afwDisplay.Display(frame=8).mtv(fitted, title="Residuals")

        fitted *= fitted
        fitted /= var

        if display:
            afwDisplay.Display(frame=9).mtv(fitted, title="Chi2")

        return fneg, negPsfSum, fpos, posPsfSum, fitted

    def testPsfDipoleFit(self, scaling=100.):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc, scaling=scaling)
        source = self.measureDipole(s, exposure)
        # Recreate the simultaneous joint Psf fit in python
        fp = source.getFootprint()
        peaks = fp.getPeaks()
        speaks = [(p.getPeakValue(), p) for p in peaks]
        speaks.sort()
        dpeaks = [speaks[0][1], speaks[-1][1]]

        negCenter = geom.Point2D(dpeaks[0].getFx(), dpeaks[0].getFy())
        posCenter = geom.Point2D(dpeaks[1].getFx(), dpeaks[1].getFy())

        fneg, negPsfSum, fpos, posPsfSum, residIm = self._makeModel(exposure, psf, fp, negCenter, posCenter)

        # Should be close to the same as the inputs; as fracOffset
        # gets smaller this will be worse.  This works for scaling =
        # 100.
        self.assertAlmostEqual(1e-2*scaling, -1e-2*fneg, 2)
        self.assertAlmostEqual(1e-2*scaling, 1e-2*fpos, 2)

        # Now compare the LeastSquares results fitted here to the C++
        # implementation: Since total flux is returned, and this is of
        # order 1e4 for this default test, scale back down so that
        # assertAlmostEqual behaves reasonably (the comparison to 2
        # places means to 0.01).  Also note that PsfDipoleFlux returns
        # the total flux, while here we are just fitting for the
        # scaling of the Psf.  Therefore the comparison is
        # fneg*negPsfSum to flux.dipole.psf.neg.
        self.assertAlmostEqual(1e-4*fneg*negPsfSum,
                               1e-4*source.get("ip_diffim_PsfDipoleFlux_neg_instFlux"),
                               2)
        self.assertAlmostEqual(1e-4*fpos*posPsfSum,
                               1e-4*source.get("ip_diffim_PsfDipoleFlux_pos_instFlux"),
                               2)

        self.assertGreater(source.get("ip_diffim_PsfDipoleFlux_pos_instFluxErr"), 0.0)
        self.assertGreater(source.get("ip_diffim_PsfDipoleFlux_neg_instFluxErr"), 0.0)
        self.assertFalse(source.get("ip_diffim_PsfDipoleFlux_neg_flag"))
        self.assertFalse(source.get("ip_diffim_PsfDipoleFlux_pos_flag"))

        self.assertAlmostEqual(source.get("ip_diffim_PsfDipoleFlux_centroid_x"), 50.0, 1)
        self.assertAlmostEqual(source.get("ip_diffim_PsfDipoleFlux_centroid_y"), 50.0, 1)
        self.assertAlmostEqual(source.get("ip_diffim_PsfDipoleFlux_neg_centroid_x"), negCenter[0], 1)
        self.assertAlmostEqual(source.get("ip_diffim_PsfDipoleFlux_neg_centroid_y"), negCenter[1], 1)
        self.assertAlmostEqual(source.get("ip_diffim_PsfDipoleFlux_pos_centroid_x"), posCenter[0], 1)
        self.assertAlmostEqual(source.get("ip_diffim_PsfDipoleFlux_pos_centroid_y"), posCenter[1], 1)
        self.assertFalse(source.get("ip_diffim_PsfDipoleFlux_neg_flag"))
        self.assertFalse(source.get("ip_diffim_PsfDipoleFlux_pos_flag"))

        self.assertGreater(source.get("ip_diffim_PsfDipoleFlux_chi2dof"), 0.0)

    def measureDipole(self, s, exp):
        msConfig = ipDiffim.DipoleMeasurementConfig()
        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("centroid_x", type=float)
        schema.addField("centroid_y", type=float)
        schema.addField("centroid_flag", type='Flag')
        task = ipDiffim.DipoleMeasurementTask(schema, config=msConfig)
        measCat = afwTable.SourceCatalog(schema)
        measCat.defineCentroid("centroid")
        source = measCat.addNew()
        source.set("centroid_x", self.xc)
        source.set("centroid_y", self.yc)
        source.setFootprint(s.getFootprint())
        # Then run the default SFM task.  Results not checked
        task.run(measCat, exp)
        return measCat[0]

    def testDipoleAnalysis(self):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        source = self.measureDipole(s, exposure)
        dpAnalysis = ipDiffim.DipoleAnalysis()
        dpAnalysis(source)

    def testDipoleDeblender(self):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        source = self.measureDipole(s, exposure)
        dpDeblender = ipDiffim.DipoleDeblender()
        dpDeblender(source, exposure)


class DipoleMeasurementTaskTest(lsst.utils.tests.TestCase):
    """A test case for the DipoleMeasurementTask.  Essentially just
    test the classification flag since the invididual algorithms are
    tested above"""

    def setUp(self):
        np.random.seed(666)
        self.config = ipDiffim.DipoleMeasurementConfig()

    def tearDown(self):
        del self.config

    def testMeasure(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = ipDiffim.DipoleMeasurementTask(schema, config=self.config)
        table = afwTable.SourceTable.make(schema)
        sources = afwTable.SourceCatalog(table)
        source = sources.addNew()
        # make fake image
        psf, psfSum, exposure, s = createDipole(100, 100, 50, 50)

        # set it in source with the appropriate schema
        source.setFootprint(s.getFootprint())
        task.run(sources, exposure)
        self.assertEqual(source.get("ip_diffim_ClassificationDipole_value"), 1.0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
