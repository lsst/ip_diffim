import lsst.ip.diffim as ipDiffim
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
import numpy as num
import lsst.afw.image.testUtils as imTestUtils

klist = ipDiffim.makeAlardLuptonBasisList(15, 3, [2, 4, 8], [4, 3, 2])
frame = 1
for kernel in klist:
    kImageOut = afwImage.ImageD(kernel.getDimensions())
    kSum      = kernel.computeImage(kImageOut, False)
    ds9.mtv(kImageOut, frame=frame)
    frame += 1

kim1 = afwImage.ImageD(klist[0].getDimensions())
kim2 = afwImage.ImageD(klist[0].getDimensions())

for k1 in range(0, len(klist)):
    klist[k1].computeImage(kim1, False)
    print k1, num.sum(num.ravel(imTestUtils.arrayFromImage(kim1)))
    
for k1 in range(0, len(klist)):
    klist[k1].computeImage(kim1, False)
    arr1 = imTestUtils.arrayFromImage(kim1)
    
    for k2 in range(k1, len(klist)):
        klist[k2].computeImage(kim2, False)
        arr2 = imTestUtils.arrayFromImage(kim2)

        #print k1, k2, num.sum(arr1*arr2)
        
        
        
