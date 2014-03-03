#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import unittest
import numpy as np
import lsst.utils.tests as tests
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
#display = True
try:
    display
except:
    display = False
np.random.seed(666)
sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def createDipole(w, h, xc, yc, scaling = 100.0, fracOffset = 1.2):
    # Make random noise image: set image plane to normal distribution
    image = afwImage.MaskedImageF(w,h)
    image.set(0)
    array = image.getImage().getArray()
    array[:,:] = np.random.randn(w,h)
    # Set variance to 1.0
    var   = image.getVariance()
    var.set(1.0)

    if display:
        ds9.mtv(image, frame=1, title="Original image")
        ds9.mtv(image.getVariance(), frame=2, title="Original variance")

    # Create Psf for dipole creation and measurement
    psfSize = 17
    psf = measAlg.DoubleGaussianPsf(psfSize, psfSize, 2.0, 3.5, 0.1)
    psfFwhmPix = sigma2fwhm * psf.computeShape().getDeterminantRadius()
    psfim = psf.computeImage().convertF()
    psfim *= scaling / psf.computePeak()
    psfw, psfh = psfim.getDimensions()
    psfSum = np.sum(psfim.getArray())
    
    # Create the dipole, offset by fracOffset of the Psf FWHM (pixels)
    offset = fracOffset * psfFwhmPix // 2
    array  = image.getImage().getArray()
    xp, yp = xc - psfw//2 + offset, yc - psfh//2 + offset
    array[yp:yp+psfh, xp:xp+psfw] += psfim.getArray()

    xn, yn = xc - psfw//2 - offset, yc - psfh//2 - offset
    array[yn:yn+psfh, xn:xn+psfw] -= psfim.getArray()

    if display:
        ds9.mtv(image, frame=3, title="With dipole")
    
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
        ds9.mtv(image, frame=4, title="Detection plane")
    
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
    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class DipoleAlgorithmTest(unittest.TestCase):
    """ A test case for dipole algorithms"""
    def setUp(self):
        self.w, self.h = 100, 100 # size of image
        self.xc, self.yc = 50, 50 # location of center of dipole

    def tearDown(self):
        pass

    def testNaiveDipoleCentroid(self):
        alg = ipDiffim.NaiveDipoleCentroidControl()

        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        schema  = afwTable.SourceTable.makeMinimalSchema()
        msb     = measAlg.MeasureSourcesBuilder()\
                   .addAlgorithm(alg)
        ms      = msb.build(schema)
        table   = afwTable.SourceTable.make(schema)
        source  = table.makeRecord()
        
        source.setFootprint(s.getFootprint())
        ms.apply(source, exposure, afwGeom.Point2D(self.xc, self.yc))
        for key in (".pos", ".pos.err", ".pos.flags", ".neg", ".neg.err", ".neg.flags"):
            try:
                source.get(alg.name+key)
            except:
                self.fail()

    def testNaiveDipoleFluxControl(self):
        alg = ipDiffim.NaiveDipoleFluxControl()

        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        schema  = afwTable.SourceTable.makeMinimalSchema()
        msb     = measAlg.MeasureSourcesBuilder()\
                   .addAlgorithm(alg)
        ms      = msb.build(schema)
        table   = afwTable.SourceTable.make(schema)
        source  = table.makeRecord()
        
        source.setFootprint(s.getFootprint())
        ms.apply(source, exposure, afwGeom.Point2D(self.xc, self.yc))
        for key in (".pos", ".pos.err", ".pos.flags", ".neg", ".neg.err", 
                    ".neg.flags", ".npos", ".nneg"):
            try:
                source.get(alg.name+key)
            except:
                self.fail()

    def testPsfDipoleFluxControl(self):
        alg = ipDiffim.PsfDipoleFluxControl()

        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        schema  = afwTable.SourceTable.makeMinimalSchema()
        msb     = measAlg.MeasureSourcesBuilder()\
                   .addAlgorithm(alg)
        ms      = msb.build(schema)
        table   = afwTable.SourceTable.make(schema)
        source  = table.makeRecord()
        
        source.setFootprint(s.getFootprint())
        ms.apply(source, exposure, afwGeom.Point2D(self.xc, self.yc))
        for key in (".pos", ".pos.err", ".pos.flags", ".neg", ".neg.err", ".neg.flags"):
            try:
                source.get(alg.name+key)
            except:
                self.fail()

    def testAll(self):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        source = self.measureDipole(s, exposure)

    def _makeModel(self, exposure, psf, fp, negCenter, posCenter):

        negPsf = psf.computeImage(negCenter).convertF()
        posPsf = psf.computeImage(posCenter).convertF()
        negPeak = psf.computePeak(negCenter)
        posPeak = psf.computePeak(posCenter)
        negPsf /= negPeak
        posPsf /= posPeak

        model    = afwImage.ImageF(fp.getBBox())
        negModel = afwImage.ImageF(fp.getBBox())
        posModel = afwImage.ImageF(fp.getBBox())

        # The center of the Psf should be at negCenter, posCenter
        negPsfBBox = negPsf.getBBox(afwImage.PARENT)
        posPsfBBox = posPsf.getBBox(afwImage.PARENT)
        modelBBox  = model.getBBox(afwImage.PARENT)

        # Portion of the negative Psf that overlaps the montage
        negOverlapBBox = afwGeom.Box2I(negPsfBBox)
        negOverlapBBox.clip(modelBBox)
        self.assertFalse(negOverlapBBox.isEmpty())

        # Portion of the positivePsf that overlaps the montage
        posOverlapBBox = afwGeom.Box2I(posPsfBBox)
        posOverlapBBox.clip(modelBBox)
        self.assertFalse(posOverlapBBox.isEmpty())

        negPsfSubim    = type(negPsf)(negPsf, negOverlapBBox, afwImage.PARENT)
        modelSubim     = type(model)(model, negOverlapBBox, afwImage.PARENT)
        negModelSubim  = type(negModel)(negModel, negOverlapBBox, afwImage.PARENT)
        modelSubim    += negPsfSubim  # just for debugging
        negModelSubim += negPsfSubim  # for fitting

        posPsfSubim    = type(posPsf)(posPsf, posOverlapBBox, afwImage.PARENT)
        modelSubim     = type(model)(model, posOverlapBBox, afwImage.PARENT)
        posModelSubim  = type(posModel)(posModel, posOverlapBBox, afwImage.PARENT)
        modelSubim    += posPsfSubim
        posModelSubim += posPsfSubim

        data = afwImage.ImageF(exposure.getMaskedImage().getImage(), fp.getBBox())
        var = afwImage.ImageF(exposure.getMaskedImage().getVariance(), fp.getBBox())
        matrixNorm = 1. / np.sqrt(np.median(var.getArray()))

        if display:
            ds9.mtv(model, frame=5, title="Unfitted model")
            ds9.mtv(data, frame=6, title="Data")

        posPsfSum = np.sum(posPsf.getArray())
        negPsfSum = np.sum(negPsf.getArray())

        M = np.array((np.ravel(negModel.getArray()), np.ravel(posModel.getArray()))).T.astype(np.float64)
        B = np.array((np.ravel(data.getArray()))).astype(np.float64)
        M *= matrixNorm
        B *= matrixNorm

        # Numpy solution
        fneg0, fpos0 = np.linalg.lstsq(M, B)[0]

        # Afw solution
        lsq = afwMath.LeastSquares.fromDesignMatrix(M, B, afwMath.LeastSquares.DIRECT_SVD)
        fneg, fpos = lsq.getSolution()

        # Should be exaxtly the same as each other
        self.assertAlmostEqual(1e-2*fneg0,  1e-2*fneg)
        self.assertAlmostEqual(1e-2*fpos0,  1e-2*fpos)
        
        # Recreate model
        fitted  = afwImage.ImageF(fp.getBBox())
        negFit  = type(negPsf)(negPsf, negOverlapBBox, afwImage.PARENT, True)
        negFit *= float(fneg)
        posFit  = type(posPsf)(posPsf, posOverlapBBox, afwImage.PARENT, True)
        posFit *= float(fpos)

        fitSubim  = type(fitted)(fitted, negOverlapBBox, afwImage.PARENT)
        fitSubim += negFit
        fitSubim  = type(fitted)(fitted, posOverlapBBox, afwImage.PARENT)
        fitSubim += posFit
        if display:
            ds9.mtv(fitted, frame=7, title="Fitted model")

        fitted   -= data

        if display:
            ds9.mtv(fitted, frame=8, title="Residuals")

        fitted   *= fitted
        fitted   /= var

        if display:
            ds9.mtv(fitted, frame=9, title="Chi2")

        return fneg, negPsfSum, fpos, posPsfSum, fitted
        
    def testPsfDipoleFit(self, scaling=100.):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc, scaling=scaling)
        source = self.measureDipole(s, exposure)

        # Recreate the simultaneous joint Psf fit in python
        fp     = source.getFootprint()
        peaks  = fp.getPeaks()
        speaks = [(p.getPeakValue(), p) for p in peaks]
        speaks.sort() 
        dpeaks = [speaks[0][1], speaks[-1][1]]

        negCenter = afwGeom.Point2D(dpeaks[0].getFx(), dpeaks[0].getFy())
        posCenter = afwGeom.Point2D(dpeaks[1].getFx(), dpeaks[1].getFy())

        fneg, negPsfSum, fpos, posPsfSum, residIm = self._makeModel(exposure, psf, fp, negCenter, posCenter)

        # Should be close to the same as the inputs; as fracOffset
        # gets smaller this will be worse.  This works for scaling =
        # 100.
        self.assertAlmostEqual(1e-2*scaling, -1e-2*fneg, 2)
        self.assertAlmostEqual(1e-2*scaling,  1e-2*fpos, 2)

        # Now compare the LeastSquares results fitted here to the C++
        # implementation: Since total flux is returned, and this is of
        # order 1e4 for this default test, scale back down so that
        # assertAlmostEqual behaves reasonably (the comparison to 2
        # places means to 0.01).  Also note that PsfDipoleFlux returns
        # the total flux, while here we are just fitting for the
        # scaling of the Psf.  Therefore the comparison is
        # fneg*negPsfSum to flux.dipole.psf.neg.
        self.assertAlmostEqual(1e-4*fneg*negPsfSum, 1e-4*source.get("flux.dipole.psf.neg"), 2)
        self.assertAlmostEqual(1e-4*fpos*posPsfSum, 1e-4*source.get("flux.dipole.psf.pos"), 2)
        self.assertTrue(source.get("flux.dipole.psf.pos.err") > 0.0)
        self.assertTrue(source.get("flux.dipole.psf.neg.err") > 0.0)
        self.assertEqual(source.get("flux.dipole.psf.neg.flags"), False)
        self.assertEqual(source.get("flux.dipole.psf.pos.flags"), False)

        self.assertAlmostEqual(source.get("flux.dipole.psf.centroid")[0], 50.0, 1)
        self.assertAlmostEqual(source.get("flux.dipole.psf.centroid")[1], 50.0, 1)
        self.assertAlmostEqual(source.get("flux.dipole.psf.neg.centroid")[0], negCenter[0], 1)
        self.assertAlmostEqual(source.get("flux.dipole.psf.neg.centroid")[1], negCenter[1], 1)
        self.assertAlmostEqual(source.get("flux.dipole.psf.pos.centroid")[0], posCenter[0], 1)
        self.assertAlmostEqual(source.get("flux.dipole.psf.pos.centroid")[1], posCenter[1], 1)
        self.assertEqual(source.get("flux.dipole.psf.centroid.flags"), False)
        self.assertEqual(source.get("flux.dipole.psf.neg.centroid.flags"), False)
        self.assertEqual(source.get("flux.dipole.psf.pos.centroid.flags"), False)

        self.assertAlmostEqual(source.get("flux.dipole.psf.chi2dof"), 1.0, 2)

        

    def measureDipole(self, s, exp):
        schema  = afwTable.SourceTable.makeMinimalSchema()
        msb     = measAlg.MeasureSourcesBuilder()\
                   .addAlgorithm(ipDiffim.NaiveDipoleCentroidControl())\
                   .addAlgorithm(ipDiffim.NaiveDipoleFluxControl())\
                   .addAlgorithm(ipDiffim.PsfDipoleFluxControl())
        ms      = msb.build(schema)
        table   = afwTable.SourceTable.make(schema)
        source  = table.makeRecord()
        source.setFootprint(s.getFootprint())
        ms.apply(source, exp, afwGeom.Point2D(self.xc, self.yc))
        return source 

    def testDipoleAnalysis(self):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        source = self.measureDipole(s, exposure)
        dpAnalysis = ipDiffim.DipoleAnalysis()
        results = dpAnalysis(source)

    def testDipoleDeblender(self):
        psf, psfSum, exposure, s = createDipole(self.w, self.h, self.xc, self.yc)
        source = self.measureDipole(s, exposure)
        dpDeblender = ipDiffim.DipoleDeblender()
        deblendSource = dpDeblender(source, exposure)
        

class DipoleMeasurementTaskTest(unittest.TestCase):
    """A test case for the DipoleMeasurementTask.  Essentially just
    test the classification flag since the invididual algorithms are
    tested above"""
    def setUp(self):
        self.config = ipDiffim.DipoleMeasurementConfig()
        self.config.algorithms.names.add("centroid.dipole.naive")
        self.config.algorithms.names.add("flux.dipole.naive")
        self.config.algorithms.names.add("flux.dipole.psf")

    def tearDown(self):
        del self.config

    def testMeasure(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        dipoleFlag = ipDiffim.DipoleMeasurementTask._ClassificationFlag
        schema.addField(dipoleFlag, "F", "probability of being a dipole")
        task = ipDiffim.DipoleMeasurementTask(schema, config=self.config)
        table = afwTable.SourceTable.make(schema)
        sources = afwTable.SourceCatalog(table)
        source = sources.addNew()
    
        # make fake image
        psf, psfSum, exposure, s = createDipole(100, 100, 50, 50)

        # set it in source with the appropriate schema
        source.setFootprint(s.getFootprint())

        task.run(exposure, sources)
        self.assertEqual(source.get(dipoleFlag), 1.0)

def suite():
    """Returns a suite containing all the test cases in this module."""

    tests.init()

    suites = []
    suites += unittest.makeSuite(DipoleAlgorithmTest)
    suites += unittest.makeSuite(DipoleMeasurementTaskTest)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

