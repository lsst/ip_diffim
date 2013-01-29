import numpy
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.ip.diffim as ipDiffim

w,h = 100,100
xc,yc=50,50

# Make random noise image
image = afwImage.MaskedImageF(w,h)
image.set(0)
array = image.getImage().getArray()
array[:,:] = numpy.random.randn(w,h)
var   = image.getVariance()
var.set(1.0)
ds9.mtv(image, frame=1)
ds9.mtv(image.getVariance(), frame=2)
# Create Psf for dipole creation and measurement
psf = afwDet.createPsf("DoubleGaussian", 11, 11, 1.0, 2.5, 0.1)
psfim = psf.computeImage(afwGeom.Point2D(0., 0.)).convertF()
psfw, psfh = psfim.getDimensions()

# For the dipole
array = image.getImage().getArray()
x0, y0 = xc-psfw+psfw//3, yc-psfh+psfh//3
array[y0:y0+psfh, x0:x0+psfw] += psfim.getArray() * 100 
x0, y0 = xc-psfw//3, yc-psfh//3
array[y0:y0+psfh, x0:x0+psfw] -= psfim.getArray() * 100
ds9.mtv(image, frame=3)

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
ds9.mtv(image, frame=4)

# Merge them together
print "# Before merge, nSrc =", len(results.sources)
fpSet = results.fpSets.positive
fpSet.merge(results.fpSets.negative, 0, 0, False)
sources = afwTable.SourceCatalog(table)
fpSet.makeSources(sources)
print "# After merge, nSrc =", len(sources)
s = sources[0]
print "# And nPeaks =", len(s.getFootprint().getPeaks())

# Measure dipole at known location
schema  = afwTable.SourceTable.makeMinimalSchema()
control = ipDiffim.NaiveDipoleCentroidControl()
msb     = measAlg.MeasureSourcesBuilder().addAlgorithm(control)
ms      = msb.build(schema)
table   = afwTable.SourceTable.make(schema)
source  = table.makeRecord()

# ms.apply(source, exp, afwGeom.Point2D(xc, yc))
# Hmm throws Assertion `px != 0' failed.
# Looks like this is because there is no getFootprint() yet...

source.setFootprint(s.getFootprint())
ms.apply(source, exp, afwGeom.Point2D(xc, yc))

# Hmm, this also fails
# print source[control.name]
# KeyError: "Field 'centroid.dipole.naive' not found in Schema."
# I think thats because the fields are 'centroid.dipole.naive.pos' and '.neg'

for key in schema.getNames():
    print key, source.get(key)


#control = ipDiffim.NaiveDipoleFluxControl()
#msb     = msb.addAlgorithm(control)

