import lsst.afw.image as afwImage
import lsst.afw.math  as afwMath
import lsst.ip.diffim as ipDiffim

cellSize  = 10
fullSize  = 100
stampSize = 10
coord     = 7.

tmi = afwImage.MaskedImageF(stampSize, stampSize)
smi = afwImage.MaskedImageF(stampSize, stampSize)

bbox    = afwImage.BBox(afwImage.PointI(0,0), fullSize, fullSize)
cellSet = afwMath.SpatialCellSet(bbox, cellSize, cellSize)
cand    = ipDiffim.makeKernelCandidate(coord, coord, tmi, smi)
cellSet.insertCandidate(cand)
