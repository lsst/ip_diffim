import lsst.afw.image as afwImage
import lsst.afw.math  as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.policy as pexPolicy

cellSize  = 10
fullSize  = 100
stampSize = 10
coord     = 7.

tmi = afwImage.MaskedImageF(stampSize, stampSize)
smi = afwImage.MaskedImageF(stampSize, stampSize)

bbox    = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(fullSize, fullSize))
cellSet = afwMath.SpatialCellSet(bbox, cellSize, cellSize)
policy  = pexPolicy.Policy()
policy.set("candidateCoreRadius", 3)
cand    = ipDiffim.makeKernelCandidate(coord, coord, tmi, smi, policy)
print cand.getStatus()
cellSet.insertCandidate(cand)
