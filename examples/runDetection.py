import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import os, sys

if len(sys.argv) == 1:
    defDataDir = eups.productDir("afwdata")
    defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                   "v5-e0-c011-a00.sci")
    templateImage  = afwImage.MaskedImageF(defTemplatePath)
else:
    templateImage  = afwImage.MaskedImageF(sys.argv[1])
    
ds9.mtv(templateImage, frame=0)

bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
bctrl.setNxSample(max(2, int(templateImage.getWidth()/256) + 1))
bctrl.setNySample(max(2, int(templateImage.getHeight()/256) + 1))
bctrl.sctrl.setNumSigmaClip(3)
bctrl.sctrl.setNumIter(3)
im      = templateImage.getImage()
backobj = afwMath.makeBackground(im, bctrl)
im     -= backobj.getImageF()

detSet         = afwDetection.makeDetectionSet(templateImage, afwDetection.createThreshold(4, "stdev"))
footprints     = detSet.getFootprints()

for fp in footprints:
    bboxes = afwDetection.footprintToBBoxList(fp)
    for bbox in bboxes:

        x0, y0, x1, y1 = bbox.getX0(), bbox.getY0(), bbox.getX1(), bbox.getY1()
        
        x0 -= 0.5; y0 -= 0.5
        x1 += 0.5; y1 += 0.5
        
        ds9.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype=ds9.RED)

bitmask = templateImage.getMask().getPlaneBitMask("DETECTED")
afwDetection.setMaskFromFootprintList(templateImage.getMask(), footprints, bitmask)
ds9.mtv(templateImage, frame=1)

