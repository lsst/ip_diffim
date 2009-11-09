import os, sys
import eups
import time
import lsst.afw.image.imageLib as afwImage
import lsst.afw.math.mathLib as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLogging

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

subBackground = True
display = True

Verbosity = 4
pexLogging.Trace_setVerbosity("lsst.ip.diffim", Verbosity)

defDataDir   = eups.productDir("afwdata") 
imageProcDir = eups.productDir("ip_diffim")

if len(sys.argv) == 1:
    defTemplatePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2_tmpl")
    defSciencePath  = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2")
elif len(sys.argv) == 3:
    defTemplatePath = sys.argv[1]
    defSciencePath  = sys.argv[2]
else:
    sys.exit(1)
    
defPolicyPath   = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")
defOutputPath   = "diffImage"

templateImage = afwImage.ExposureF(defTemplatePath)
scienceImage  = afwImage.ExposureF(defSciencePath)
policy        = ipDiffim.generateDefaultPolicy(defPolicyPath)

if defSciencePath.find('CFHT') == -1:
    templateImage = ipDiffim.warpTemplateExposure(templateImage,
                                                  scienceImage,
                                                  policy)
    policy.set("detThreshold", 5.)
    policy.set("sizeCellX", 128)
    policy.set("sizeCellY", 128)
    policy.set("kernelBasisSet", "alard-lupton")
    policy.set("useRegularization", False)

if subBackground:
    # NOTE - you need to subtract off background from the image
    # you run detection on.  Here it is the template.
    algorithm   = policy.get("backgroundPolicy.algorithm")
    binsize     = policy.get("backgroundPolicy.binsize")
    undersample = policy.get("backgroundPolicy.undersample")
    bctrl       = afwMath.BackgroundControl(algorithm)
    bctrl.setNxSample(templateImage.getWidth()//binsize + 1)
    bctrl.setNySample(templateImage.getHeight()//binsize + 1)
    bctrl.setUndersampleStyle(undersample)
    
    image   = templateImage.getMaskedImage().getImage() 
    backobj = afwMath.makeBackground(image, bctrl)
    image  -= backobj.getImageF()
    
    image   = scienceImage.getMaskedImage().getImage() 
    backobj = afwMath.makeBackground(image, bctrl)
    image  -= backobj.getImageF()
    
    del image; del backobj
    #
    #

spatialKernel, spatialBg, kernelCellSet = ipDiffim.createPsfMatchingKernel(templateImage.getMaskedImage(),
                                                                           scienceImage.getMaskedImage(),
                                                                           policy)

# Lets see what we got
if display:
    mos = displayUtils.Mosaic()

    # Inputs
    frame = 0
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False): # False = include bad candidates
            cand  = ipDiffim.cast_KernelCandidateF(cand)
            rchi2 = cand.getChi2()
                
            try:
                im = cand.getImage()
                if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    statStr = "Good"
                elif cand.getStatus() == afwMath.SpatialCellCandidate.BAD:
                    statStr = "Bad"
                else:
                    statStr = "Unkn"
                mos.append(im, "#%d: %.1f (%s)" % (cand.getId(), rchi2, statStr))
            except Exception, e:
                pass
    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)

    # Bases
    frame = 1
    mos.reset()
    basisList = spatialKernel.getKernelList()
    for idx in range(len(basisList)):
        kernel = basisList[idx]
        im   = afwImage.ImageD(spatialKernel.getDimensions())
        ksum = kernel.computeImage(im, False)
        mos.append(im, "K%d" % (idx))
    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)
        

    # Spatial model
    frame = 2
    mos.reset()
    width = templateImage.getWidth()
    height = templateImage.getHeight()
    stamps = []; stampInfo = []
    for x in (0, width//2, width):
        for y in (0, height//2, height):
            im   = afwImage.ImageD(spatialKernel.getDimensions())
            ksum = spatialKernel.computeImage(im,
                                              False,
                                              afwImage.indexToPosition(x),
                                              afwImage.indexToPosition(y))
            mos.append(im, "x=%d y=%d kSum=%.2f" % (x, y, ksum))

    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)
            

    # Background
    frame = 3
    backgroundIm = afwImage.ImageF(templateImage.getWidth(),
                                   templateImage.getHeight(), 0)
    ipDiffim.addToImage(backgroundIm, spatialBg)
    ds9.mtv(backgroundIm, frame=frame)

    # Diffim!
    frame = 4
    diffIm = ipDiffim.convolveAndSubtract(templateImage.getMaskedImage(),
                                          scienceImage.getMaskedImage(),
                                          spatialKernel,
                                          spatialBg)
    ds9.mtv(diffIm, frame=frame)


# examples/runSpatialModel.py $AFWDATA_DIR/DC3a-Sim/sci/v5-e0/v5-e0-c011-a00.sci $AFWDATA_DIR/DC3a-Sim/sci/v26-e0/v26-e0-c011-a00.sci
