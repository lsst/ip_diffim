import os
import re
import sys

import lsst.afw.image as afwImage

TESTING = False

# of the form /share/lsst9/DC2root/RUN0002/ipd/input/705310/01/705310p_01_img.fits
inputExposureRoot = sys.argv[1]
inputExposurePath = os.path.abspath(inputExposureRoot)
inputExposureDirs = inputExposurePath.split('/')

inputExposureFile = inputExposureDirs[-1]
inputExposureAmp  = inputExposureDirs[-2]
inputExposureID   = inputExposureDirs[-3]
inputExposureBase = '/'.join(inputExposureDirs[:-3])

if not TESTING:
    inputMaskedImage = afwImage.MaskedImageF()
    inputMaskedImage.readFits(inputExposureRoot)
    inputWCS = afwImage.Wcs(inputMaskedImage.getImage().getMetaData())
    inputExposure = afwImage.ExposureF(inputMaskedImage, inputWCS)

    nRowSubexposures = 4 # int(sys.argv[2]) # 4
    nColSubexposures = 2 # int(sys.argv[3]) # 2
    nRowMaskedImage = inputMaskedImage.getRows() # 4644
    nColMaskedImage = inputMaskedImage.getCols() # 2112
    
else:
    nRowSubexposures = 4
    nColSubexposures = 2
    nRowMaskedImage  = 4644
    nColMaskedImage  = 2112

nRowPix = int(nRowMaskedImage / nRowSubexposures)
nColPix = int(nColMaskedImage / nColSubexposures)

extn = 0
for row in range(nRowSubexposures):
    for col in range(nColSubexposures):

        outputExposureBase = os.path.join(inputExposureBase, '%s%d' % (inputExposureID, extn))
        if not os.path.isdir(outputExposureBase):
            cmd = 'mkdir %s' % (outputExposureBase)
            print '#', cmd
            if not TESTING:
                os.system(cmd)

        outputExposureAmp = os.path.join(outputExposureBase, inputExposureAmp)
        if not os.path.isdir(outputExposureAmp):
            cmd = 'mkdir %s' % (outputExposureAmp)
            print '#', cmd
            if not TESTING:
                os.system(cmd)
            

        outputExposureFile = re.sub(inputExposureID, '%s%d' % (inputExposureID, extn), inputExposureFile)
        outputExposureFile = os.path.join(outputExposureAmp, outputExposureFile)
        print '# Writing', outputExposureFile
        
        bbox = afwImage.BBox2i(col * nColPix,
                         row * nRowPix,
                         nColPix,
                         nRowPix)

        if not TESTING:
            outputExposure = inputExposure.getSubExposure(bbox)
            #outputExposure = inputMaskedImage.getSubImage(bbox)
            outputExposure.writeFits(outputExposureFile)
            
        extn += 1
