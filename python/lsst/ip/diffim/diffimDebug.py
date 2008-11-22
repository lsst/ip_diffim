import lsst.ip.diffim.diffimPlot as ipDiffimPlot
import lsst.ip.diffim.diffimTools as ipDiffimTools
import numpy
from lsst.pex.logging import Trace
import lsst.daf.base as dafBase

def plotDiffImQuality1(diffIm, kernel, template, image, label, outfile=None):
    data     = ipDiffimTools.imageToMatrix(diffIm.getImage())
    variance = ipDiffimTools.imageToMatrix(diffIm.getVariance())
    mask     = ipDiffimTools.imageToMatrix(diffIm.getMask())
    idx      = numpy.where(mask == 0)
    sigma    = numpy.ravel( data[idx] / numpy.sqrt(variance[idx]) )
    Trace('lsst.ip.diffim', 5,
          '%s : Python diffim residuals = %.2f +/- %.2f sigma' % (label, sigma.mean(), sigma.std()))

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
    Trace('lsst.ip.diffim', 5,
          'Kernel %d : Python diffim residuals = %.2f +/- %.2f sigma' % (id, cSigma.mean(), cSigma.std()))
   
    dData     = ipDiffimTools.imageToMatrix(dDiffIm.getImage())
    dVariance = ipDiffimTools.imageToMatrix(dDiffIm.getVariance())
    dMask     = ipDiffimTools.imageToMatrix(dDiffIm.getMask())
    dIdx      = numpy.where(dMask == 0)
    dSigma    = numpy.ravel( dData[dIdx] / numpy.sqrt(dVariance[dIdx]) )
    Trace('lsst.ip.diffim', 5,
          'Kernel %d : Python diffim residuals = %.2f +/- %.2f sigma' % (id, dSigma.mean(), dSigma.std()))
    
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


def writeDiffImages(id,
                    tStamp, iStamp,
                    cDifi, cDiffIm, cKernel, 
                    dDifi, dDiffIm, dKernel):
    
    iStamp.writeFits('iFoot_%d' % (id))
    tStamp.writeFits('tFoot_%d' % (id))
    
    ckp,cks = cKernel.computeNewImage(False)
    cmd = ckp.getMetaData()
    cmd.addProperty(dafBase.DataProperty('CONV', 'Template'))
    cmd.addProperty(dafBase.DataProperty('MSIG', cDifi.getSingleStats().getResidualMean()))
    cmd.addProperty(dafBase.DataProperty('VSIG', cDifi.getSingleStats().getResidualStd()))
    cmd.addProperty(dafBase.DataProperty('KSUM', cks))
    ckp.setMetadata(cmd)
    ckp.writeFits('cKernel_%d.fits' % (id))
    
    dkp,dks = dKernel.computeNewImage(False)
    dmd = dkp.getMetaData()
    dmd.addProperty(dafBase.DataProperty('CONV', 'Template'))
    dmd.addProperty(dafBase.DataProperty('MSIG', dDifi.getSingleStats().getResidualMean()))
    dmd.addProperty(dafBase.DataProperty('VSIG', dDifi.getSingleStats().getResidualStd()))
    dmd.addProperty(dafBase.DataProperty('KSUM', dks))
    dkp.setMetadata(dmd)
    dkp.writeFits('dKernel_%d.fits' % (id))
    
    cDiffIm.writeFits('cDiff_%d' % (id))
    dDiffIm.writeFits('dDiff_%d' % (id))
    
    #cSigma = ipDiffimTools.imageToMatrix(convDiffIm2.getImage()) / numpy.sqrt(ipDiffimTools.imageToMatrix(convDiffIm2.getVariance()))
    #cSigma = ipDiffimTools.matrixToImage(cSigma)
    #cSigma.writeFits('cSig_%d.fits' % (footprintID))
    
    #dcSigma = ipDiffimTools.imageToMatrix(deconvDiffIm2.getImage()) / numpy.sqrt(ipDiffimTools.imageToMatrix(deconvDiffIm2.getVariance()))
    #dcSigma = ipDiffimTools.matrixToImage(dcSigma)
    #dcSigma.writeFits('dcSig_%d.fits' % (footprintID))
    
