import lsst.afw.math  as afwMath
import lsst.afw.image as afwImage
import numpy
from lsst.pex.logging import Trace
import lsst.daf.base as dafBase
from   lsst.pex.logging import Trace

# relative imports, since these are in __init__.py
import diffimLib 
import diffimDebug
import diffimPlot 

import pdb

def plotDiffImQuality1(difi, diffim, kernel, label, outfile=None):
    template = difi.getImageToConvolvePtr().get()
    image    = difi.getImageToNotConvolvePtr().get()

    data     = imageToMatrix(diffim.getImage())
    variance = imageToMatrix(diffim.getVariance())
    mask     = imageToMatrix(diffim.getMask())
    idx      = numpy.where(mask == 0)
    sigma    = numpy.ravel( data[idx] / numpy.sqrt(variance[idx]) )

    info     = (imageToMatrix(template.getImage()),
                imageToMatrix(image.getImage()),
                data / numpy.sqrt(variance),
                imageToMatrix(kernel.computeNewImage(False)[0]),
                sigma)
    ipDiffimPlot.plotSigmaHistograms( (info,), title=label, outfile=outfile )

  
    

# Debugging plots etc.
def plotDiffImQuality2(id, iteration,
                       cDiffIm, cKernel, cTemplate, cImage,
                       dDiffIm, dKernel, dTemplate, dImage):
    # Need to do this to make the pixel histogram
    cData     = imageToMatrix(cDiffIm.getImage())
    cVariance = imageToMatrix(cDiffIm.getVariance())
    cMask     = imageToMatrix(cDiffIm.getMask())
    cIdx      = numpy.where(cMask == 0)
    cSigma    = numpy.ravel( cData[cIdx] / numpy.sqrt(cVariance[cIdx]) )
   
    dData     = imageToMatrix(dDiffIm.getImage())
    dVariance = imageToMatrix(dDiffIm.getVariance())
    dMask     = imageToMatrix(dDiffIm.getMask())
    dIdx      = numpy.where(dMask == 0)
    dSigma    = numpy.ravel( dData[dIdx] / numpy.sqrt(dVariance[dIdx]) )
    
    cInfo     = (imageToMatrix(cTemplate.getImage()),
                 imageToMatrix(cImage.getImage()),
                 cData / numpy.sqrt(cVariance),
                 imageToMatrix(cKernel.computeNewImage(False)[0]),
                 cSigma)
    
    dInfo     = (imageToMatrix(dTemplate.getImage()),
                 imageToMatrix(dImage.getImage()),
                 dData / numpy.sqrt(dVariance),
                 imageToMatrix(dKernel.computeNewImage(False)[0]),
                 dSigma)
    ipDiffimPlot.plotSigmaHistograms( (cInfo, dInfo), title='Kernel %d' % (id), outfile='Kernel_%d_%d.ps' % (id, iteration) )


def writeDiffImages(prefix, id, kModel):

    if not kModel.isBuilt():
        Trace('lsst.ip.diffim.diffimDebug.writeDiffImages', 5, 'Building model for debugging')
        kModel.buildModel()

#    if not kModel.getStatus():
#        # I need to be careful here; since I am debuggin I want to
#        # look at everything except for when there was a kernel
#        # exception.
#        return
        
    kModel.getMiToConvolvePtr().writeFits('tFoot_%s_%s' % (prefix, id))
    kModel.getMiToNotConvolvePtr().writeFits('iFoot_%s_%s' % (prefix, id))

    cki = afwImage.ImageD(kModel.getKernelPtr().getDimensions())
    cks = kModel.getKernelPtr().computeImage(cki, False)
    cmd = cki.getMetadata()
    if cmd == None:
        cmd = dafBase.PropertySet()

    cmd.setString('CONV', 'Template')
    cmd.setFloat('KCOL', kModel.getColc())
    cmd.setFloat('KROW', kModel.getRowc())
#    cmd.setFloat('MSIG', kModel.getStats().getMean())
#    cmd.setFloat('VSIG', kModel.getStats().getVariance())
    cmd.setFloat('BG',   kModel.getBackground())
    cmd.setFloat('KSUM', cks)
    cmd.setBool('KQUALITY', kModel.getStatus())
    cki.setMetadata(cmd)
    cki.writeFits('kernel_%s_%s.fits' % (prefix, id))
    
    diffIm = diffimLib.convolveAndSubtract(kModel.getMiToConvolvePtr(),
                                           kModel.getMiToNotConvolvePtr(),
                                           kModel.getKernelPtr(),
                                           kModel.getBackground())
    
    diffIm.writeFits('diff_%s_%s' % (prefix, id))



#######
# Temporary functions until we formalize this in the build system somewhere
#######

def imageToMatrix(im, dtype=float):
    arr = numpy.zeros([im.getCols(), im.getRows()], dtype=dtype)
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            arr[col, row] = im.getVal(col, row)
    return arr

def imageToVector(im, dtype=float):
    arr = numpy.zeros([im.getCols()*im.getRows()], dtype=dtype)
    n   = 0
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            arr[n] = im.getVal(col, row)
            n += 1
    return arr

def matrixToImage(arr):
    im = afwImage.ImageF(arr.shape[0], arr.shape[1])
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[col, row])
    return im

def vectorToImage(arr, nCols, nRows):
    assert len(arr) == nCols * nRows
    im = afwImage.ImageF(nCols, nRows)
    n  = 0
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[n])
            n += 1
    return im

def matrixToKernelPtr(arr):
    im = afwImage.ImageD(arr.shape[0], arr.shape[1])
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[col, row])
    return afwMath.KernelPtr( afwMath.FixedKernel(im) ) 

def vectorToKernelPtr(arr, nCols, nRows):
    # need imageD for FixedKernels
    assert len(arr) == nCols * nRows
    im = afwImage.ImageD(nCols, nRows)
    n  = 0
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[n])
            n += 1
    return afwMath.KernelPtr( afwMath.FixedKernel(im) ) 

def vectorPairToVectors(vectorPair):
    kernelVector = afwMath.vectorD()
    kernelErrorVector = afwMath.vectorD()
    for i in range(vectorPair.size()-1):
        kernelVector.push_back(vectorPair[i][0])
        kernelErrorVector.push_back(vectorPair[i][1])

    background = vectorPair.back()[0]
    backgroundError = vectorPair.back()[1]
    return kernelVector, kernelErrorVector, background, backgroundError
