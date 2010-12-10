import os
import sys
import lsst.afw.image.imageLib as afwImage
import lsst.afw.math.mathLib as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as pexLogging

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

subBackground = True
display = True
fwhm = 6.8
warp = True

verbosity = 5
pexLogging.Trace_setVerbosity("lsst.ip.diffim", verbosity)

defDataDir   = eups.productDir("afwdata") 
imageProcDir = eups.productDir("ip_diffim")

if len(sys.argv) == 1:
    defSciencePath  = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2")
    defTemplatePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2_tmpl")
elif len(sys.argv) == 3:
    defSciencePath  = sys.argv[1]
    defTemplatePath = sys.argv[2]
else:
    sys.exit(1)
    
defOutputPath   = "diffImage"

templateImage = afwImage.ExposureF(defTemplatePath)
scienceImage  = afwImage.ExposureF(defSciencePath)
policy        = ipDiffim.createDefaultPolicy()

if warp:
    templateImage = ipDiffim.warpTemplateExposure(templateImage,
                                                  scienceImage,
                                                  policy.getPolicy("warpingPolicy"))
    #policy.set("detThreshold", 5.)
    #policy.set("sizeCellX", 128)
    #policy.set("sizeCellY", 128)
    #policy.set("kernelBasisSet", "alard-lupton")
    #policy.set("useRegularization", False)

if subBackground:
    diffimTools.backgroundSubtract(policy.getPolicy("afwBackgroundPolicy"),
                                   [templateImage.getMaskedImage(),
                                    scienceImage.getMaskedImage()])

frame  = 0
ds9.mtv(templateImage, frame=frame)
frame += 1
ds9.mtv(scienceImage, frame=frame)

spatialKernel, spatialBg, kernelCellSet = ipDiffim.makePsfMatchingKernel(templateImage.getMaskedImage(),
                                                                         scienceImage.getMaskedImage(),
                                                                         policy)

import pdb
pdb.set_trace()

# Lets see what we got
if display:
    mos = displayUtils.Mosaic()

    # Inputs
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
    frame += 1
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)

    # KernelCandidates
    frame += 1
    diffimTools.displayCandiateResults(kernelCellSet, frame)

    # Bases
    mos.reset()
    basisList = spatialKernel.getKernelList()
    for idx in range(len(basisList)):
        kernel = basisList[idx]
        im   = afwImage.ImageD(spatialKernel.getDimensions())
        ksum = kernel.computeImage(im, False)
        mos.append(im, "K%d" % (idx))
    mosaic = mos.makeMosaic()
    frame += 1
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)
        

    # Spatial model
    mos.reset()
    width = templateImage.getWidth()
    height = templateImage.getHeight()
    stamps = []
    stampInfo = []
    for x in (0, width//2, width):
        for y in (0, height//2, height):
            im   = afwImage.ImageD(spatialKernel.getDimensions())
            ksum = spatialKernel.computeImage(im,
                                              False,
                                              afwImage.indexToPosition(x),
                                              afwImage.indexToPosition(y))
            mos.append(im, "x=%d y=%d kSum=%.2f" % (x, y, ksum))

    mosaic = mos.makeMosaic()
    frame += 1
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)
            

    # Background
    backgroundIm = afwImage.ImageF(templateImage.getWidth(),
                                   templateImage.getHeight(), 0)
    backgroundIm += spatialBg
    frame += 1
    ds9.mtv(backgroundIm, frame=frame)

    # Diffim!
    diffIm = ipDiffim.convolveAndSubtract(templateImage.getMaskedImage(),
                                          scienceImage.getMaskedImage(),
                                          spatialKernel,
                                          spatialBg)
    frame += 1
    ds9.mtv(diffIm, frame=frame)


# examples/runSpatialModel.py $AFWDATA_DIR/DC3a-Sim/sci/v5-e0/v5-e0-c011-a00.sci
# ... $AFWDATA_DIR/DC3a-Sim/sci/v26-e0/v26-e0-c011-a00.sci
