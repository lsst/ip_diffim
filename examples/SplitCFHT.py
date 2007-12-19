import lsst.fw.Core.fwLib as fw
import sys

inputExposureRoot = sys.argv[1]

inputMaskedImage = fw.MaskedImageF()
inputMaskedImage.readFits(inputExposureRoot)
inputWCS = fw.WCS(inputMaskedImage.getImage().getMetaData())
inputExposure = fw.ExposureF(inputMaskedImage, inputWCS)

nRowSubexposures = int(sys.argv[2]) # 4
nColSubexposures = int(sys.argv[3]) # 2

nRowMaskedImage = inputMaskedImage.getRows() # 4644
nColMaskedImage = inputMaskedImage.getCols() # 2112

print nRowMaskedImage, nColMaskedImage

nRowPix = int(nRowMaskedImage / nRowSubexposures)
nColPix = int(nColMaskedImage / nColSubexposures)

extn = 0
for row in range(nRowSubexposures):
    for col in range(nColSubexposures):
        bbox = fw.BBox2i(col * nColPix,
                         row * nRowPix,
                         nColPix,
                         nRowPix)


        outputExposure = inputExposure.getSubExposure(bbox)
        outputImage = inputMaskedImage.getSubImage(bbox)
        outputFilename = '%s%d' % (inputExposureRoot, extn)
        print extn, ':', bbox.min().x(), bbox.max().x(), bbox.min().y(), bbox.max().y(), outputFilename
        #outputExposure.writeFits(outputFilename)
        extn += 1
