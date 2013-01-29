import numpy
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg

w,h = 100,100
image = afwImage.MaskedImageF(w,h)
image.set(0)

xc,yc=50,50
array = image.getImage().getArray()
array[:,:] = numpy.random.randn(w,h)
ds9.mtv(image, frame=1)


# For the dipole
psf = afwDet.createPsf("DoubleGaussian", 11, 11, 1.0, 2.5, 0.1)
psfim = psf.computeImage(afwGeom.Point2D(0., 0.)).convertF()
psfim *= 100
psfw, psfh = psfim.getDimensions()

array = image.getImage().getArray()
array[yc-psfw//2:yc-psfw//2+psfw, xc-psfw//2:xc-psfw//2+psfw] += psfim.getArray()
array[yc-psfw//4:yc-psfw//4+psfh, xc-psfw//4:xc-psfw//4+psfw] -= psfim.getArray()
      
ds9.mtv(image, frame=2)

exp = afwImage.makeExposure(image)
exp.setPsf(psf)
config = measAlg.SourceDetectionConfig()
config.thresholdPolarity = "both"
config.reEstimateBackground = False
schema = afwTable.SourceTable.makeMinimalSchema()
task = measAlg.SourceDetectionTask(schema, config=config)
table = afwTable.SourceTable.make(schema)
results = task.makeSourceCatalog(table, exp)
ds9.mtv(image, frame=3)

print len(results.sources)
fpSet = results.fpSets.positive
fpSet.merge(results.fpSets.negative, 0, 0, False)
sources = afwTable.SourceCatalog(table)
fpSet.makeSources(sources)
print len(sources)
s = sources[0]
print len(s.getFootprint().getPeaks())





