#!/usr/bin/env python
import os, pdb, sys, re

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.policy as pexPolicy
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as logging
import numpy as num
import lsst.afw.display.ds9 as ds9

display = True

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

# Evaluates DC3a processing

class DC3aTestCase():

    def __init__(self):
        self.setUp()
    
    def setUp(self):
        self.policy      = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')
        self.basisList   = ipDiffim.generateDeltaFunctionKernelSet(self.kCols, self.kRows)
        self.minPix      = self.policy.getInt('fpNpixMin')
        self.maxPix      = self.policy.getInt('fpNpixMax')
        self.fpGrowSize  = self.kCols # int( self.policy.getDouble('fpGrowKsize') * self.kCols )

        # difference imaging functor
        self.kFunctor    = ipDiffim.PsfMatchingFunctorF(self.basisList)

        # image statistics
        self.dStats      = ipDiffim.ImageStatisticsF()

        # background subtraction
        self.bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
        self.bctrl.sctrl.setNumSigmaClip(3)
        self.bctrl.sctrl.setNumIter(3)

    def reInitialize(self, sciPath, tmplPath, diffPath, detSig=10):
        self.scienceImage     = afwImage.ExposureF(sciPath)
        self.differenceImage  = afwImage.ExposureF(diffPath)
        self.differenceImage.getMaskedImage().setXY0( self.scienceImage.getMaskedImage().getXY0() )

        self.bctrl.setNxSample(max(2, int(self.scienceImage.getWidth()/256) + 1))
        self.bctrl.setNySample(max(2, int(self.scienceImage.getHeight()/256) + 1))
        im           = self.scienceImage.getMaskedImage().getImage()
        self.backobj = afwMath.makeBackground(im, self.bctrl)
        im          -= self.backobj.getImageF()

        threshold   = afwDetection.createThreshold(5, 'variance')
        self.detSet = afwDetection.makeDetectionSet(self.scienceImage.getMaskedImage(), threshold)

        #print detSet.getFootprints().size()
        #self.footprints = detSet.getFootprints()
        #print self.footprints.size()
        print '#', sciPath, self.detSet.getFootprints().size()

        # Ok, because of ticket #835 we can't trust the tmpl/ images
        # on disk.  And reading the monolithic template is memory
        # intensive.  So we will only search the difference image
        # here, and hope that Ray can produce a run with Trace=5 to
        # get the "before" diffim stats.

        return
        
        self.templateImage    = afwImage.ExposureF(tmplPath)

        # Remap the template to the image; replace self.templateImage with warped image
        wKernel = afwMath.makeWarpingKernel('lanczos4')
        self.remappedImage = self.templateImage.Factory(
            self.scienceImage.getWidth(), 
            self.scienceImage.getHeight(),
            self.scienceImage.getWcs())
        self.remappedImage.getMaskedImage().setXY0( self.scienceImage.getMaskedImage().getXY0() )

        # Hack due to ticket #835
        filter           = self.scienceImage.getMetadata().get('FILTER').strip()
        datasetId        = self.scienceImage.getMetadata().get('datasetId').strip()

        #templatePath     = os.path.join('/lsst/images/repository/template', datasetId, '25', filter, 'T0004_%s_25_%s' % (datasetId, filter))
        templatePath     = os.path.join('/lsst/becker/lsst_devel/DMS/afwdata_trunk/templates', 'T0004_%s_25_%s' % (datasetId, filter))
        
        templateExposure = afwImage.ExposureF(templatePath)
        afwMath.warpExposure(self.remappedImage, 
                             templateExposure, 
                             wKernel)
        self.templateImage = self.remappedImage

        # footprints
        self.footprints = ipDiffim.getCollectionOfFootprintsForPsfMatching(
            self.templateImage.getMaskedImage(),
            self.scienceImage.getMaskedImage(),
            self.policy)
        print '#', len(self.footprints)

        if display:
            ds9.mtv(self.scienceImage, frame=0)
            ds9.mtv(self.templateImage, frame=1)
            ds9.mtv(self.differenceImage, frame=2)
            sys.exit(1)
        
    def tearDown(self):
        del self.policy

    def applyFunctor(self, footprint):
        # check for any masked pixels
        sBits = ipDiffim.FindSetBitsU(self.scienceImage.getMaskedImage().getMask())
        sBits.apply( footprint )
        if sBits.getBits() > 0:
            return None, None

        sBits = ipDiffim.FindSetBitsU(self.templateImage.getMaskedImage().getMask())
        sBits.apply( footprint )
        if sBits.getBits() > 0:
            return None, None
        
        try:
            smi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  footprint.getBBox())
            tmi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), footprint.getBBox())
        except:
            return None

        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        # accepts : imageToConvolve, imageToNotConvolve
        try:
            self.kFunctor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
        except:
            return None, None
        kernel    = self.kFunctor.getKernel()
        diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.reset()
        self.dStats.apply( diffIm2 )

        return self.dStats.getMean(), self.dStats.getRms()

    def getStats2(self, footprint, count):
        sBits = ipDiffim.FindSetBitsU(self.differenceImage.getMaskedImage().getMask())

        # Will throw if off image
        try:
            sBits.apply( footprint )
        except:
            #print '# Fail 1'
            return None
        
        if sBits.getBits() > 0:
            #print '# Fail 2'
            return None

        bbox  = footprint.getBBox()
        bbox.shift(-self.differenceImage.getMaskedImage().getX0(),
                   -self.differenceImage.getMaskedImage().getY0())
        try:
            smi = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(), bbox)
            dmi = afwImage.MaskedImageF(self.differenceImage.getMaskedImage(), bbox)
        except:
            #print '# Fail 3'
            return None

        smi.writeFits('SMI'+str(count))
        dmi.writeFits('DMI'+str(count))

        srcFlux  = afwMath.makeStatistics(smi.getImage(),    afwMath.SUM).getValue()
        diffMean = afwMath.makeStatistics(dmi.getImage(),    afwMath.MEAN).getValue()
        diffVar  = afwMath.makeStatistics(dmi.getImage(),    afwMath.VARIANCE).getValue()
        varMean  = afwMath.makeStatistics(dmi.getVariance(), afwMath.MEAN).getValue()

        return srcFlux, diffMean, diffVar, varMean

    def getStats(self, footprint):
        sBits = ipDiffim.FindSetBitsU(self.differenceImage.getMaskedImage().getMask())

        # Will throw if off image
        try:
            sBits.apply( footprint )
        except:
            #print '# Fail 1'
            return None
        
        if sBits.getBits() > 0:
            #print '# Fail 2'
            return None

        bbox  = footprint.getBBox()
        bbox.shift(-self.differenceImage.getMaskedImage().getX0(),
                   -self.differenceImage.getMaskedImage().getY0())
        try:
            smi = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(), bbox)
            dmi = afwImage.MaskedImageF(self.differenceImage.getMaskedImage(), bbox)
        except:
            #print '# Fail 3'
            return None

        flux = afwMath.makeStatistics(smi.getImage(), afwMath.SUM).getValue()
        self.dStats.reset()
        self.dStats.apply( dmi )

        return self.dStats.getMean(), self.dStats.getRms(), flux
    
def run(ntodo):
    """Run the tests"""

    count = 0

    myInfo = []

    dc3a = DC3aTestCase()
    
    #root = '/lsst/becker/lsst_devel/DMS/afwdata_trunk/DC3a-CFHT'
    root = '/lsst/DC3root/rlp1233/IPSD/output/'
    # use sci/ as the base directory
    sci  = os.path.join(root, 'sci')
    tmpl = os.path.join(root, 'tmpl')
    diff = os.path.join(root, 'diff')

    count = 0
    for dir1 in os.listdir(sci):
        path1 = os.path.join(sci, dir1)
        
        for file in os.listdir(path1):
            sci_img = os.path.join(path1, file)
            
            if not sci_img.endswith('_img.fits'):
                continue
    
            sci_msk = re.sub('_img', '_msk', sci_img)
            sci_var = re.sub('_img', '_var', sci_img)
            if not os.path.exists(sci_msk):
                continue
            if not os.path.exists(sci_var):
                continue
    
            tmpl_root = dir1.split('-')[0]
            tmpl_e    = dir1.split('-')[1]
            tmpl_img  = os.path.join(tmpl, tmpl_root, re.sub('.sci_img', '.tmpl_img', re.sub('-'+tmpl_e, '', file)))
            tmpl_msk  = os.path.join(tmpl, tmpl_root, re.sub('.sci_img', '.tmpl_msk', re.sub('-'+tmpl_e, '', file)))
            tmpl_var  = os.path.join(tmpl, tmpl_root, re.sub('.sci_img', '.tmpl_var', re.sub('-'+tmpl_e, '', file)))
            if not os.path.exists(tmpl_img):
                continue
            if not os.path.exists(tmpl_msk):
                continue
            if not os.path.exists(tmpl_var):
                continue
    
            diff_img  = os.path.join(diff, dir1, re.sub('.sci_img', '.diff_img', file))
            diff_msk  = os.path.join(diff, dir1, re.sub('.sci_img', '.diff_msk', file))
            diff_var  = os.path.join(diff, dir1, re.sub('.sci_img', '.diff_var', file))
            if not os.path.exists(diff_img):
                continue
            if not os.path.exists(diff_msk):
                continue
            if not os.path.exists(diff_var):
                continue
    
            # do it!
            dc3a.reInitialize( re.sub('_img.fits', '', sci_img),
                               re.sub('_img.fits', '', tmpl_img),
                               re.sub('_img.fits', '', diff_img) )

            footprints = dc3a.detSet.getFootprints()
            for fpId in range(footprints.size()):
                fp = footprints[fpId]
                if (fp.getNpix() < dc3a.minPix) or (fp.getNpix() > dc3a.maxPix):
                    continue
                # Use this for the statistics
                fpGrow = afwDetection.growFootprint(fp, dc3a.fpGrowSize, False)

                fpX = int( 0.5*(fpGrow.getBBox().getX0()+fpGrow.getBBox().getX1()) )
                fpY = int( 0.5*(fpGrow.getBBox().getY0()+fpGrow.getBBox().getY1()) )
                bgX = fpX - dc3a.scienceImage.getMaskedImage().getX0()
                bgY = fpY - dc3a.scienceImage.getMaskedImage().getY0()

                results = dc3a.getStats2(fpGrow, count)
                count += 1
                
                if results == None:
                    continue
                srcFlux, diffMean, diffVar, varMean = results
                print fpId, fpX, fpY, srcFlux, diffMean, diffVar, varMean

                #results = dc3a.getStats(fpGrow)
                #if results == None:
                #    continue
                #mean, rms, flux = results
                #bg    = dc3a.backobj.getPixel(bgX, bgY)
                #bgTot = bg * fpGrow.getNpix()
                #print fpId, fpX, fpY, mean, rms, flux, bgTot
                #myInfo.append( (mean, rms, flux, bgTot) )

            count += 1
            if count == ntodo:
                return myInfo

    return myInfo


def myhist(means, rms, output):
    pylab.figure()
    
    sp1 = pylab.subplot(121)
    sp1.hist(means, bins=50)
    #sp1.set_title('%.3f +/- %.3f' % (means.mean(), means.std()))
    sp1.set_xlabel('Mean residual in sigma', fontsize=12, weight='bold')
    sp1.set_ylabel('N objects', fontsize=12, weight='bold')
    pylab.setp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize=10)
    
    sp2 = pylab.subplot(122)
    sp2.hist(rms,   bins=50)
    #sp2.set_title('%.3f +/- %.3f' % (rms.mean(), rms.std()))
    sp2.set_xlabel('Std of residuals in sigma', fontsize=12, weight='bold')
    pylab.setp(sp2.get_xticklabels()+sp2.get_yticklabels(), fontsize=10)
    
    pylab.savefig(output)

    
if __name__ == "__main__":
    ntodo = int(sys.argv[1])
    results = run(ntodo)
    sys.exit(1)
    
    import pylab
    means = num.array( [ x[0] for x in results ] )
    rms   = num.array( [ x[1] for x in results ] )
    flux  = num.array( [ x[2] for x in results ] )
    bg    = num.array( [ x[3] for x in results ] )

    frat = flux / bg
    med  = num.median(frat)
    idx1 = num.where(frat < med)
    idx2 = num.where(frat > med)

    print '# All', means.mean(), means.std(), rms.mean(), rms.std()
    print '# Faint', means[idx1].mean(), means[idx1].std(), rms[idx1].mean(), rms[idx1].std()
    print '# Bright', means[idx2].mean(), means[idx2].std(), rms[idx2].mean(), rms[idx2].std()
    
    myhist(means, rms, 'plotAll.png')
    myhist(means[idx1], rms[idx1], 'plotFaint.png')
    myhist(means[idx2], rms[idx2], 'plotBright.png')
