import eups
import sys, os, optparse
import lsst.daf.base  as dafBase
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.afw.detection as detection
from lsst.pex.logging import Log
from lsst.pex.logging import Trace
from lsst.pex.policy  import Policy
import pdb
import time
import numpy as num
import pylab

def testBasics(templateMaskedImage, scienceMaskedImage, policy, kSize=19, gScale=1):
    policy.set('kernelCols', kSize)
    policy.set('kernelRows',  kSize)

    kImage     = afwImage.ImageD(kSize, kSize)
    kBasisList = ipDiffim.generateDeltaFunctionKernelSet(kSize, kSize)
    kFunctor   = ipDiffim.PsfMatchingFunctorF(kBasisList)

    # the Eigen matrix with be this_size x this_size
    nKernelPars = 1 + kBasisList.size()
    print 'Matrices are', nKernelPars, 'x', nKernelPars

    # a footprint of size 1, grown by kSize, will end up being
    # (2*kSize+1)  x  (2*kSize+1)
    fp     = detection.Footprint(afwImage.BBox(afwImage.PointI(484, 525), 1, 1))
    fpGrow = detection.growFootprint(fp, int(gScale*kSize))
    fpBBox = fpGrow.getBBox()

    tSubImage = afwImage.MaskedImageF(templateMaskedImage, fpBBox)
    iSubImage = afwImage.MaskedImageF(scienceMaskedImage, fpBBox)

    # the good pixels post-convolution in kFunctor.apply()
    startCol = kBasisList[0].getCtrX()
    startRow = kBasisList[0].getCtrY()
    endCol   = tSubImage.getWidth()  - (kBasisList[0].getWidth()  - kBasisList[0].getCtrX()) + 1
    endRow   = tSubImage.getHeight() - (kBasisList[0].getHeight() - kBasisList[0].getCtrY()) + 1
    nConstraints = (startCol-endCol) * (startRow-endRow)
    print 'Constraints number', nConstraints

    # if a single-pixel footprint is grown by kSize-1 pixels, we will
    # have N-1 constraints, where N = the matrices we need to invert.
    #
    # so the minimum number of pixels a single-pixel footprint can be
    # grown is kSize

    model = ipDiffim.SpatialModelKernelF(fpGrow,
                                         tSubImage,
                                         iSubImage,
                                         kFunctor,
                                         policy,
                                         False)
    t1 = time.time()
    model.buildModel()
    t2 = time.time()
    return model.getKernelSum(), model.getBg(), model.getBgErr(), t2-t1

def checkSizes(templateMaskedImage, scienceMaskedImage, policy, xctr, yctr, kSize=19, gScale=1):
    policy.set('kernelCols', kSize)
    policy.set('kernelRows',  kSize)

    kImage     = afwImage.ImageD(kSize, kSize)
    kBasisList = ipDiffim.generateDeltaFunctionKernelSet(kSize, kSize)
    kFunctor   = ipDiffim.PsfMatchingFunctorF(kBasisList)

    # the Eigen matrix with be this_size x this_size
    nKernelPars = 1 + kBasisList.size()
    print 'Matrices are', nKernelPars, 'x', nKernelPars
    
        
    # a footprint of size 1, grown by kSize, will end up being
    # (2*kSize+1)  x  (2*kSize+1)
    fp     = detection.Footprint(afwImage.BBox(afwImage.PointI(xctr, yctr), 1, 1))
    fpGrow = detection.growFootprint(fp, int(gScale*kSize))
    fpBBox = fpGrow.getBBox()
    
    tSubImage = afwImage.MaskedImageF(templateMaskedImage, fpBBox)
    iSubImage = afwImage.MaskedImageF(scienceMaskedImage, fpBBox)
    
    # the good pixels post-convolution in kFunctor.apply()
    startCol = kBasisList[0].getCtrX()
    startRow = kBasisList[0].getCtrY()
    endCol   = tSubImage.getWidth()  - (kBasisList[0].getWidth()  - kBasisList[0].getCtrX()) + 1
    endRow   = tSubImage.getHeight() - (kBasisList[0].getHeight() - kBasisList[0].getCtrY()) + 1
    nConstraints = (startCol-endCol) * (startRow-endRow)
    print 'Constraints number', nConstraints
    
    # if a single-pixel footprint is grown by kSize-1 pixels, we will
    # have N-1 constraints, where N = the matrices we need to invert.
    #
    # so the minimum number of pixels a single-pixel footprint can be
    # grown is kSize
    
    model = ipDiffim.SpatialModelKernelF(fpGrow,
                                         tSubImage,
                                         iSubImage,
                                         kFunctor,
                                         policy,
                                         False)
    t1 = time.time()
    model.buildModel()
    t2 = time.time()
    return model.getKernelSum(), model.getBg(), model.getBgErr(), t2-t1

    
def main():
    imageProcDir = eups.productDir('ip_diffim')
    if imageProcDir == None:
        print 'Error: could not set up ip_diffim'
        sys.exit(1)
    defPolicyPath   = os.path.join(imageProcDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
    defVerbosity    = 0
    
    usage = ""
    parser = optparse.OptionParser(usage)
    parser.add_option('-p', '--policy', default=defPolicyPath, help='policy file')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
        help='verbosity of diagnostic trace messages; 1 for just warnings, more for more information')
    
    (options, args) = parser.parse_args()

    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue

    sciencePath  = args[0]
    templatePath = args[1]
    outputPath   = args[2]
    policyPath   = options.policy

    print 'Science image: ', sciencePath
    print 'Template image:', templatePath
    print 'Output image:  ', outputPath
    print 'Policy file:   ', policyPath

    templateMaskedImage = afwImage.MaskedImageF(templatePath)
    scienceMaskedImage  = afwImage.MaskedImageF(sciencePath)
    policy              = Policy.createPolicy(policyPath)
    
    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)

    fpList = ipDiffim.getCollectionOfFootprintsForPsfMatching(
        templateMaskedImage,
        scienceMaskedImage,
        policy)

    for fpID, fpPtr in enumerate(fpList):
        xctr = int(0.5 * (fpPtr.getBBox().getX0() + fpPtr.getBBox().getX1()))
        yctr = int(0.5 * (fpPtr.getBBox().getY0() + fpPtr.getBBox().getY1()))
        
        arrays = []
        scales = num.arange(1, 3.1, 0.1)
        for i in scales:
            arrays.append( checkSizes(templateMaskedImage, scienceMaskedImage, policy, xctr, yctr, gScale=i) )

        ksums  = num.array( [ x[0] for x in arrays ] )
        bgs    = num.array( [ x[1] for x in arrays ] )
        dbgs   = num.array( [ x[2] for x in arrays ] )
        ts     = num.array( [ x[3] for x in arrays ] )
        
        pylab.figure()
        
        sp1 = pylab.subplot(221)
        sp1.plot(scales, bgs/dbgs, 'r-')
        sp1.set_ylabel('Background / Bg Error', fontsize=10)
        sp1.set_xlabel('Grow Footprint by X Kernel Sizes', fontsize=10)
        sp1.set_title('%d %d' % (xctr, yctr))
                         
        sp2 = pylab.subplot(222)
        sp2.plot(scales, dbgs, 'r-')
        sp2.set_ylabel('Background Error', fontsize=10)
        sp2.set_xlabel('Grow Footprint', fontsize=10)
        
        sp3 = pylab.subplot(223)
        sp3.plot(scales, ksums, 'r-')
        sp3.set_ylabel('Kernel Sums', fontsize=10)
        sp3.set_xlabel('Grow Footprint', fontsize=10)
        
        sp4 = pylab.subplot(224)
        sp4.plot(scales, ts, 'r-')
        sp4.set_ylabel('Processing time (s)', fontsize=10)
        sp4.set_xlabel('Grow Footprint', fontsize=10)
        
        pylab.setp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize=8)
        pylab.setp(sp2.get_xticklabels()+sp2.get_yticklabels(), fontsize=8)
        pylab.setp(sp3.get_xticklabels()+sp3.get_yticklabels(), fontsize=8)
        pylab.setp(sp4.get_xticklabels()+sp4.get_yticklabels(), fontsize=8)
        
        #    sp1.semilogy()
        sp2.semilogy()
        #    sp3.semilogy()
        #    sp4.semilogy()
        
        sp1.set_xlim( scales[0]-0.1, scales[-1]+0.1 )
        sp2.set_xlim( scales[0]-0.1, scales[-1]+0.1 )
        sp3.set_xlim( scales[0]-0.1, scales[-1]+0.1 )
        sp4.set_xlim( scales[0]-0.1, scales[-1]+0.1 )

        pylab.savefig('fp%d.png' % (fpID))
        
    pylab.show()
    
def run():
    Log.getDefaultLog()
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), 'Objects leaked:'
        print dafBase.Citizen_census(dafBase.cout, memId0)

if __name__ == '__main__':
    run()
