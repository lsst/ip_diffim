import lsst.ip.diffim.diffimPlot as ipDiffimPlot
import lsst.ip.diffim.diffimTools as ipDiffimTools
import numpy
from lsst.pex.logging import Trace
import lsst.daf.base as dafBase

def plotDiffImQuality1(difi, diffim, kernel, label, outfile=None):
    template = difi.getImageToConvolvePtr().get()
    image    = difi.getImageToNotConvolvePtr().get()

    data     = ipDiffimTools.imageToMatrix(diffim.getImage())
    variance = ipDiffimTools.imageToMatrix(diffim.getVariance())
    mask     = ipDiffimTools.imageToMatrix(diffim.getMask())
    idx      = numpy.where(mask == 0)
    sigma    = numpy.ravel( data[idx] / numpy.sqrt(variance[idx]) )

    info     = (ipDiffimTools.imageToMatrix(template.getImage()),
                ipDiffimTools.imageToMatrix(image.getImage()),
                data / numpy.sqrt(variance),
                ipDiffimTools.imageToMatrix(kernel.computeNewImage(False)[0]),
                sigma)
    ipDiffimPlot.plotSigmaHistograms( (info,), title=label, outfile=outfile )

  
    

# Debugging plots etc.
def plotDiffImQuality2(id, iteration,
                       cDiffIm, cKernel, cTemplate, cImage,
                       dDiffIm, dKernel, dTemplate, dImage):
    # Need to do this to make the pixel histogram
    cData     = ipDiffimTools.imageToMatrix(cDiffIm.getImage())
    cVariance = ipDiffimTools.imageToMatrix(cDiffIm.getVariance())
    cMask     = ipDiffimTools.imageToMatrix(cDiffIm.getMask())
    cIdx      = numpy.where(cMask == 0)
    cSigma    = numpy.ravel( cData[cIdx] / numpy.sqrt(cVariance[cIdx]) )
   
    dData     = ipDiffimTools.imageToMatrix(dDiffIm.getImage())
    dVariance = ipDiffimTools.imageToMatrix(dDiffIm.getVariance())
    dMask     = ipDiffimTools.imageToMatrix(dDiffIm.getMask())
    dIdx      = numpy.where(dMask == 0)
    dSigma    = numpy.ravel( dData[dIdx] / numpy.sqrt(dVariance[dIdx]) )
    
    cInfo     = (ipDiffimTools.imageToMatrix(cTemplate.getImage()),
                 ipDiffimTools.imageToMatrix(cImage.getImage()),
                 cData / numpy.sqrt(cVariance),
                 ipDiffimTools.imageToMatrix(cKernel.computeNewImage(False)[0]),
                 cSigma)
    
    dInfo     = (ipDiffimTools.imageToMatrix(dTemplate.getImage()),
                 ipDiffimTools.imageToMatrix(dImage.getImage()),
                 dData / numpy.sqrt(dVariance),
                 ipDiffimTools.imageToMatrix(dKernel.computeNewImage(False)[0]),
                 dSigma)
    ipDiffimPlot.plotSigmaHistograms( (cInfo, dInfo), title='Kernel %d' % (id), outfile='Kernel_%d_%d.ps' % (id, iteration) )


def writeDiffImages(prefix, id, kModel):
    
    kModel.getMiToNotConvolvePtr().writeFits('iFoot_%s_%s' % (prefix, id))
    kModel.getMiToConvolvePtr().writeFits('tFoot_%s_%s' % (prefix, id))
    
    ckp,cks = kModel.getKernelPtr().computeNewImage(False)
    cmd = ckp.getMetaData()
    cmd.addProperty(dafBase.DataProperty('CONV', 'Template'))
    cmd.addProperty(dafBase.DataProperty('KCOL', kModel.getColcNorm()))
    cmd.addProperty(dafBase.DataProperty('KROW', kModel.getRowcNorm()))
    cmd.addProperty(dafBase.DataProperty('MSIG', kModel.getStats().getResidualMean()))
    cmd.addProperty(dafBase.DataProperty('VSIG', kModel.getStats().getResidualStd()))
    cmd.addProperty(dafBase.DataProperty('BG',   kModel.getBackground()))
    cmd.addProperty(dafBase.DataProperty('KQUALITY', kModel.getQaStatus()))
    cmd.addProperty(dafBase.DataProperty('KSUM', cks))
    ckp.setMetadata(cmd)
    ckp.writeFits('kernel_%s_%s.fits' % (prefix, id))
    
    diffIm = ipDiffim.convolveAndSubtract(kModel.getMiToConvolvePtr().get(),
                                          kModel.getMiToNotConvolvePtr().get(),
                                          kModel.getKernelPtr(),
                                          kModel.getBackground())
    
    
    diffIm.writeFits('diff_%s_%s' % (prefix, id))
